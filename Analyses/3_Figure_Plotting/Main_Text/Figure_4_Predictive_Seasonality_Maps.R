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
source("Analyses/1_Covariate_Extraction_and_Collation/2_Polygon_Construction.R")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
source("Functions/Logistic_Regression_Functions.R")
source("Functions/Spatial_Prediction_Functions.R")
India_Surroundings_Extent <- extent(60, 100, 0, 40) # Creates extent covering India - use to crop raster
source("Analyses/1_Covariate_Extraction_and_Collation/1_Loading_Rasters.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data_full <- mosquito_data
mosquito_data <- mosquito_data[, 24:35]
colnames(mosquito_data) <- seq(1, 12)
clusters <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
prior <- "informative"
set.seed(58) 
india <- getData("GADM", country = "IND", level = 0)

# Loading in Probabilities and Covariate Raster Stack
temporal_cluster_probabilities <- readRDS("Outputs/Logistic_Regression_Output/Multinomial_Logistic_Regression_Predictive_Outputs.rds")
cluster_1_prob <- temporal_cluster_probabilities$cluster_1
cluster_2_prob <- temporal_cluster_probabilities$cluster_2
cluster_3_prob <- temporal_cluster_probabilities$cluster_3
cluster_4_prob <- temporal_cluster_probabilities$cluster_4
raster_stack <- stack("Datasets/Extracted_Covariates_for_Modelling/Raster_Stack.grd")

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
threshold <- 0.66

extent <- extent(raster_stack[[1]])
cluster_1_matrix <- matrix(cluster_1_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_1_raster <- raster(cluster_1_matrix)
extent(cluster_1_raster) <- extent
cluster_1_raster <- mask(cluster_1_raster, india)
cluster_1_raster <- cluster_1_raster > threshold

cluster_2_matrix <- matrix(cluster_2_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_2_raster <- raster(cluster_2_matrix)
extent(cluster_2_raster) <- extent
cluster_2_raster <- mask(cluster_2_raster, india)
cluster_2_raster <- cluster_2_raster > threshold

cluster_3_matrix <- matrix(cluster_3_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_3_raster <- raster(cluster_3_matrix)
extent(cluster_3_raster) <- extent
cluster_3_raster <- mask(cluster_3_raster, india)
cluster_3_raster <- cluster_3_raster > threshold

cluster_4_matrix <- matrix(cluster_4_prob, nrow = 891, ncol = 890, byrow = TRUE)
cluster_4_raster <- raster(cluster_4_matrix)
extent(cluster_4_raster) <- extent
cluster_4_raster <- mask(cluster_4_raster, india)
cluster_4_raster <- cluster_4_raster > threshold

cluster_1_points <- species_geos[species_geos$cluster == 1, ]
cluster_2_points <- species_geos[species_geos$cluster == 2, ]
cluster_3_points <- species_geos[species_geos$cluster == 3, ]
cluster_4_points <- species_geos[species_geos$cluster == 4, ]

# 8.5 x 7.28
pdf(file = "Figures/Figure_6/Figure_6_Reduced_Subset_Seasonal_Profile_Predictive_Maps.pdf", height = 7.28, width = 8.5)
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
dev.off()


