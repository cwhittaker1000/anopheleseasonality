#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
##                                                                                                   ##
#######################################################################################################
library(MESS); library(numbers); library(factoextra); library(rgl); library(tsne); library(zoo); 
library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); library(shinystan); 
library(ggplot2); library(reshape2); library(deSolve); library(parallel); library(matlib); library(matlab); 
library(pracma); library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(dplyr); 
library(VGAM); library(rgl); library(car); library(caret); library(mltools); library(data.table); library(glmnet);
library(locfit); library(oce)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/")
source("R_Scripts/1_Covariate_Extraction_and_Collation/2_Polygon_Construction.R")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
source("Functions/Logistic_Regression_Functions.R")
source("Functions/Spatial_Prediction_Functions.R")
India_Surroundings_Extent <- extent(60, 100, 0, 40) # Creates extent covering India - use to crop raster
source("R_Scripts/1_Covariate_Extraction_and_Collation/1_Loading_Rasters.R")
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data_full <- mosquito_data
mosquito_data <- mosquito_data[, 24:35]
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
if (prior == "informative") {
  clusters <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
} else if (prior == "uninformative") {
  clusters <- readRDS("Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Clustering.rds")
}
covariates <- readRDS("Datasets/Extracted_Covariates_for_Modelling/Extracted_Covariates.rds")

set.seed(58) 
india <- getData("GADM", country = "IND", level = 0)

#######################################################################################################
##                                                                                                   ##
##              Loading In Fitted Negative Binomial GP Results, Normalising & Removing               ##
##                                      Low Count Time Series                                        ## 
##                                                                                                   ##
#######################################################################################################
fitted_storage <- matrix(nrow = 272, ncol = 36)
timepoints_storage <- matrix(nrow = 272, ncol = 36)
for (i in 1:272) { 
  temp <- mean_realisation_extract(i, mosquito_data, prior, FALSE)
  fitted_storage[i, 1:36] <- temp$mean
  timepoints_storage[i, 1:36] <- temp$timepoints
}
normalised <- matrix(nrow = 272, ncol = 36)
for (i in 1:272) {
  temp <-  normalise_total(fitted_storage[i, ])
  normalised[i, 1:length(temp)] <- temp
}

normalised_mosquito_data <- normalised


########################################################################################################
##                                                                                                    ##
##                       Stacking All the Rasters Together For Value Extraction                       ##
##                                                                                                    ##
##           Projecting all rasters to the same resolution and stacking them together.                ##
##                                                                                                    ##
########################################################################################################
covariates <- readRDS("Datasets/Extracted_Covariates_for_Modelling/Extracted_Covariates.rds")
unique_landcover_values <- unique(covariates[, "Dominant_Landcover"])
raster_stack <- stack("Datasets/Extracted_Covariates_for_Modelling/Raster_Stack.grd")
number_of_covariates <- nlayers(raster_stack) - 1 # to ensure dominant landcover doesn't getscaled later on
raster_value_storage <- matrix(nrow = ncell(raster_stack[[1]]), ncol = number_of_covariates + 1)
colnames(raster_value_storage) <- names(raster_stack)
for (i in 1:nlayers(raster_stack)) {
  raster_value_storage[, i] <- as.vector(raster_stack[[i]])
}
for (j in 1:length(raster_value_storage[, number_of_covariates + 1])) {
  if (!(raster_value_storage[j, number_of_covariates + 1] %in% unique_landcover_values)) {
    raster_value_storage[j, number_of_covariates + 1] <- NA
  }
}
table(raster_value_storage[, number_of_covariates + 1], useNA = "ifany")

raster_value_storage <- as.data.frame(raster_value_storage)
raster_value_storage$Dominant_Landcover <- as.factor(raster_value_storage$Dominant_Landcover)
raster_value_storage <- as.matrix(one_hot(as.data.table(raster_value_storage), cols = "Dominant_Landcover"))
number_raster_values <- dim(raster_value_storage)[2] 
new_names <- names(raster_stack)[-length(names(raster_stack))]
for (i in 1:(number_raster_values - number_of_covariates)) {
  new_names <- c(new_names, paste0("LC", i))
}
colnames(raster_value_storage) <- new_names

# Even More Reduced Subset (3 + 4 + 3 + 3 + 1 + 11 landcover variables) - 25 variables
temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Mean_Temp_Driest_Quartest")
rain_red <- c("Annual_Rain", "Rain_Seasonality", "CHIRPS_Min", "Rain_Coldest_Quarter")
arid_red <- c("Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
hydro_red <- c("Water_Areas_Occurrence", "Water_Areas_Recurrence", "Flow_Accumulation")
landcover_red <- c("City_Accessibility")
reduced_subset <- c(temp_red, rain_red, arid_red, hydro_red, landcover_red)
raster_value_storage <- cbind(raster_value_storage[, reduced_subset], raster_value_storage[, 56:66])
raster_value_storage <- as.data.frame(raster_value_storage)


########################################################################################################
##                                                                                                    ##
##                          Processing STAN Output to Retrieve Coeffecients                           ##
##                                                                                                    ##
########################################################################################################
if (prior == "informative") {
  fit <- readRDS("Outputs/Logistic_Regression_Output/Informative_Prior/Reduced_Subset_STAN_Output.rds")
} else if (prior == "uninformative") {
  fit <- readRDS("Outputs/Logistic_Regression_Output/Uninformative_Prior/Reduced_Subset_STAN_Output.rds")
}

spec_coefs <- rbind(apply(fit$alpha_ann, 2, mean), apply(fit$alpha_cul, 2, mean), apply(fit$alpha_dir, 2, mean), 
                    apply(fit$alpha_fluv, 2, mean), apply(fit$alpha_min, 2, mean), apply(fit$alpha_ste, 2, mean), apply(fit$alpha_sub, 2, mean))
row.names(spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
colnames(spec_coefs) <- c("Clus_1", "Clus_2", "Clus_3", "Clus_4")

variable_names <- c(reduced_subset, "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8", "LC9", "LC10", "LC11")
number_covariates <- length(variable_names)
envt_betas <- matrix(nrow = number_covariates, ncol = 4)
for (i in 1:number_covariates) {
  envt_betas[i, ] <- apply(fit$envt_betas[ , i, ], 2, mean)
}
row.names(envt_betas) <- variable_names
colnames(envt_betas) <- c("Clus_1", "Clus_2", "Clus_3", "Clus_4")


########################################################################################################
##                                                                                                    ##
##                     Loading In and Processing Vector Predicted Distribution Maps                   ##
##                                                                                                    ##                                                                                                    ##
########################################################################################################
ann_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.ann(complex)_ensemble.gri")
ann_distribution <- ann_distribution[[2]]
ann_distribution <- projectRaster(ann_distribution, raster_stack[[1]])
ann_distribution <- as.vector(ann_distribution)/1000
ann_distribution[is.na(ann_distribution)] <- 0

cul_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.cul_ensemble.gri")
cul_distribution <- cul_distribution[[2]]
cul_distribution <- projectRaster(cul_distribution, raster_stack[[1]])
cul_distribution <- as.vector(cul_distribution)/1000
cul_distribution[is.na(cul_distribution)] <- 0

dir_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.dir(complex)_ensemble.gri")
dir_distribution <- dir_distribution[[2]]
dir_distribution <- projectRaster(dir_distribution, raster_stack[[1]])
dir_distribution <- as.vector(dir_distribution)/1000
dir_distribution[is.na(dir_distribution)] <- 0

fluv_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.flu_ensemble.gri")
fluv_distribution <- fluv_distribution[[2]]
fluv_distribution <- projectRaster(fluv_distribution, raster_stack[[1]])
fluv_distribution <- as.vector(fluv_distribution)/1000
fluv_distribution[is.na(fluv_distribution)] <- 0

min_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.min(complex)_ensemble.gri")
min_distribution <- min_distribution[[2]]
min_distribution <- projectRaster(min_distribution, raster_stack[[1]])
min_distribution <- as.vector(min_distribution)/1000
min_distribution[is.na(min_distribution)] <- 0

ste_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.ste_ensemble.gri")
ste_distribution <- ste_distribution[[2]]
ste_distribution <- projectRaster(ste_distribution, raster_stack[[1]])
ste_distribution <- as.vector(ste_distribution)/1000
ste_distribution[is.na(ste_distribution)] <- 0

sub_distribution <- raster::stack("Datasets/Vector_Spatial_Distributions/proj_Current_Ano.sub(complex)_ensemble.gri")
sub_distribution <- sub_distribution[[2]]
sub_distribution <- projectRaster(sub_distribution, raster_stack[[1]])
sub_distribution <- as.vector(sub_distribution)/1000
sub_distribution[is.na(sub_distribution)] <- 0


########################################################################################################
##                                                                                                    ##
##        Generating Model Predictions of Seasonal Profile for Each Species and Location              ##
##                                                                                                    ##                                                                                                    ##
########################################################################################################
num_raster_cells <- ncell(raster_stack[[1]])
cluster_1_prob <- vector(mode = "numeric", length = num_raster_cells)
cluster_2_prob <- vector(mode = "numeric", length = num_raster_cells)
cluster_3_prob <- vector(mode = "numeric", length = num_raster_cells)
cluster_4_prob <- vector(mode = "numeric", length = num_raster_cells)
species_coefficient_matrix <- spec_coefs
environmental_coefficients <- as.matrix(envt_betas)
for (i in 1:num_raster_cells) {
  species_raster_values <- as.vector(c(ann_distribution[i], cul_distribution[i], dir_distribution[i], fluv_distribution[i], min_distribution[i], ste_distribution[i], sub_distribution[i]))
  environmental_raster_values <- as.matrix(raster_value_storage[i, ], nrow = 1)
  all_species_preds <- generate_all_species_predictions(species_coefficient_matrix, environmental_coefficients, species_raster_values, environmental_raster_values)
  prob_presence <- generate_probability_temporal_presence(all_species_preds)
  cluster_1_prob[i] <- prob_presence[1]
  cluster_2_prob[i] <- prob_presence[2]
  cluster_3_prob[i] <- prob_presence[3]
  cluster_4_prob[i] <- prob_presence[4]
  if (i %% 100000 == 0) {
    print(i)
  }
}

temporal_cluster_probabilities <- list(cluster_1 = cluster_1_prob, cluster_2 = cluster_2_prob, cluster_3 = cluster_3_prob, cluster_4 = cluster_4_prob)
saveRDS(temporal_cluster_probabilities, "Outputs/Logistic_Regression_Output/Informative_Prior/Reduced_Subset_Predictive_Outputs.rds")
temporal_cluster_probabilities <- readRDS("Outputs/Logistic_Regression_Output/Informative_Prior/Reduced_Subset_Predictive_Outputs.rds")
cluster_1_prob <- temporal_cluster_probabilities$cluster_1
cluster_2_prob <- temporal_cluster_probabilities$cluster_2
cluster_3_prob <- temporal_cluster_probabilities$cluster_3
cluster_4_prob <- temporal_cluster_probabilities$cluster_4


#######################################################################################################
##                                                                                                   ##
##                        Plotting the Collated Time Series On a Map of India                        ## 
##                                                                                                   ##
#######################################################################################################
par(mfrow = c(1, 1))
species_geos <- data.frame(x = NA, y = NA, species = NA, cluster = NA, stringsAsFactors = FALSE)
counter <- 1
for (i in 1:272) {
  
  index <- mosquito_data_full$Location_ID[i]
  matching_location_index <- which(location_coordinates$Location_ID %in% index)
  
  if (length(matching_location_index) == 1) {
    coords <- st_coordinates(location_coordinates$geometry[matching_location_index])[, 1:2]
    if (length(coords) == 2) {
      centroid <- coords
    } else {
      centroid <- geosphere::centroid(coords)
    }
    species_geos[counter, 1] <- centroid[1]
    species_geos[counter, 2] <- centroid[2]
    species_geos[counter, 3] <- mosquito_data_full$Species[i]
    species_geos[counter, 4] <- clusters$Cluster[i]
    counter <- counter + 1
    
  } else {
    print(length(matching_location_index))
    for (j in 1:length(matching_location_index)) {
      coords <- st_coordinates(location_coordinates$geometry[matching_location_index[j]])[, 1:2]
      if (length(coords) == 2) {
        centroid <- coords
      } else {
        centroid <- geosphere::centroid(coords)
      }
      species_geos[counter, 1] <- centroid[1]
      species_geos[counter, 2] <- centroid[2]
      species_geos[counter, 3] <- mosquito_data_full$Species[i]
      species_geos[counter, 4] <- clusters$Cluster[i]
      counter <- counter + 1
    }
  }
}


########################################################################################################
##                                                                                                    ##
##                    Masking and Plotting the Output for Cluster Probabilities                       ##
##                                                                                                    ##                                                                                                    ##
########################################################################################################
extent <- extent(raster_stack[[1]])
cluster_1_matrix <- matrix(cluster_1_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_1_raster <- raster(cluster_1_matrix)
extent(cluster_1_raster) <- extent
cluster_1_raster <- mask(cluster_1_raster, india)

cluster_2_matrix <- matrix(cluster_2_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_2_raster <- raster(cluster_2_matrix)
extent(cluster_2_raster) <- extent
cluster_2_raster <- mask(cluster_2_raster, india)

cluster_3_matrix <- matrix(cluster_3_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_3_raster <- raster(cluster_3_matrix)
extent(cluster_3_raster) <- extent
cluster_3_raster <- mask(cluster_3_raster, india)

cluster_4_matrix <- matrix(cluster_4_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_4_raster <- raster(cluster_4_matrix)
extent(cluster_4_raster) <- extent
cluster_4_raster <- mask(cluster_4_raster, india)

cluster_1_points <- species_geos[species_geos$cluster == 1, ]
cluster_2_points <- species_geos[species_geos$cluster == 2, ]
cluster_3_points <- species_geos[species_geos$cluster == 3, ]
cluster_4_points <- species_geos[species_geos$cluster == 4, ]

# 8.5 x 7.28
palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C", "#E71D36", "#52AD9C", "#7761B5", "#3F220F", "#D6D84F", "#363537"))
timepoints <- seq(1, 36)
par(mfrow = c(2, 2), oma = c(6, 6, 6, 6), mar = c(2, 0, 0, 0))
plot(cluster_1_raster, xlim = c(65, 100), ylim = c(6, 37), las = 1, legend = FALSE, xaxt = "n")
points(cluster_1_points$x, cluster_1_points$y, cex = 1, pch = 21, bg = palette()[1])
plotInset(85, 8, 97.5, 17, expr = plot(timepoints, apply(normalised_mosquito_data[clusters$Cluster == 1, ], 2, mean) * 100, type = "l", lwd = 7, col = palette()[1], xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 7)), mar=c(0, 0, 0, 0))
text(-9, 96, "A", font = 2, cex = 2.5)

plot(cluster_2_raster, xlim = c(65, 100), ylim = c(6, 37), las = 1, legend = FALSE, axes = FALSE)
points(cluster_2_points$x, cluster_2_points$y, cex = 1, pch = 21, bg = palette()[2])
plotInset(85, 8, 97.5, 17, expr = plot(timepoints, apply(normalised_mosquito_data[clusters$Cluster == 2, ], 2, mean) * 100, type = "l", lwd = 7, col = palette()[2], xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 7)), mar=c(0, 0, 0, 0))
text(-9, 96, "B", font = 2, cex = 2.5)

plot(cluster_3_raster, xlim = c(65, 100), ylim = c(6, 37), las = 1, legend = FALSE)
points(cluster_3_points$x, cluster_3_points$y, cex = 1, pch = 21, bg = palette()[3])
plotInset(85, 8, 97.5, 17, expr = plot(timepoints, apply(normalised_mosquito_data[clusters$Cluster == 3, ], 2, mean) * 100, type = "l", lwd = 7, col = palette()[3], xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 7)), mar=c(0, 0, 0, 0))
text(-9, 96, "C", font = 2, cex = 2.5)

plot(cluster_4_raster, xlim = c(65, 100), ylim = c(6, 37), las = 1, legend = FALSE, yaxt = "n")
points(cluster_4_points$x, cluster_4_points$y, cex = 1, pch = 21, bg = palette()[4])
plotInset(85, 8, 97.5, 17, expr = plot(timepoints, apply(normalised_mosquito_data[clusters$Cluster == 4, ], 2, mean) * 100, type = "l", lwd = 7, col = palette()[4], xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 7)), mar=c(0, 0, 0, 0))
text(-9, 96, "D", font = 2, cex = 2.5)

mtext("Latitude (Degrees)", side = 1, outer = TRUE, cex = 1.25, font = 2, line = 1.2, col = "grey20")
mtext("Longitude (Degrees)", side = 2, outer = TRUE, cex = 1.25, font = 2, line = 2.7, col = "grey20")
