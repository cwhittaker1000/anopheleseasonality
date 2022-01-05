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
library(UpSetR); library(dendextend); library(here)

# Sourcing functions
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
source("Functions/Logistic_Regression_Functions.R")

# Loading in and processing mosquito data
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 

if (prior == "informative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
} else if (prior == "uninformative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Clustering.rds")
}
covariates <- readRDS("Datasets/Extracted_Covariates_for_Modelling/Extracted_Covariates.rds")

#######################################################################################################
##                                                                                                   ##
##        Loading In Environmental Covariates and Cluster Assignments for Each Time Series           ## 
##                                                                                                   ##
#######################################################################################################
# Even More Reduced Subset (3 + 4 + 3 + 3 + 1 + 11 landcover variables) - 25 variables
temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Mean_Temp_Driest_Quarter")
rain_red <- c("Annual_Rain", "Rain_Seasonality", "CHIRPS_Min", "Rain_Coldest_Quarter")
arid_red <- c("Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
hydro_red <- c("Water_Areas_Occurrence", "Water_Areas_Recurrence", "Flow_Accumulation")
landcover_red <- c("City_Accessibility", "Dominant_Landcover")
reduced_subset <- c(temp_red, rain_red, arid_red, hydro_red, landcover_red)

colnames(covariates)[10] <- "Mean_Temp_Driest_Quarter"
covariates <- covariates[order(covariates[, "location_IDs"]), ]
covariates <- covariates[, -1]
covariates <- as.data.frame(covariates)
covariates <- covariates[, reduced_subset]
covariates$Dominant_Landcover <- as.factor(covariates$Dominant_Landcover)
covariates <- as.matrix(one_hot(as.data.table(covariates), cols = "Dominant_Landcover"))
number_covariates <- dim(covariates)[2] 

number_time_series <- length(cluster_labels$Cluster)
temp <- cluster_labels$Cluster
temp <- as.data.frame(cbind(Overall_Cluster = temp, Cluster = temp))
temp$Cluster <- as.factor(temp$Cluster)
Cluster_Labels <- as.data.frame(one_hot(as.data.table(temp), cols = "Cluster"))

# Creating the Datasets
Envt_Variables <- matrix(nrow = number_time_series, ncol = number_covariates)
for (i in 1:number_time_series) {
  index <- cluster_labels$Loc_ID[i]
  Envt_Variables[i, ] <- unname(as.vector(covariates[index, ]))
}
variable_names <- c(reduced_subset[-length(reduced_subset)], "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8", "LC9", "LC10", "LC11")
colnames(Envt_Variables) <- variable_names


#######################################################################################################
##                                                                                                   ##
##          Running Logistic Regression Models Implemented Using STAN With Subsampling               ## 
##                                                                                                   ##
#######################################################################################################
Number_Clusters <- length(unique(cluster_labels$Cluster))
Number_Envt_Variables <- number_covariates
Number_Obs_Ann <- sum(species == "Annularis")
Number_Obs_Cul <- sum(species == "Culicifacies")
Number_Obs_Dir <- sum(species == "Dirus")
Number_Obs_Fluv <- sum(species == "Fluviatilis")
Number_Obs_Min <- sum(species == "Minimus")
Number_Obs_Ste <- sum(species == "Stephensi")
Number_Obs_Sub <- sum(species == "Subpictus")

options(mc.cores = 3)
STAN_GLM <- stan_model("Model_Files/Multinomial_Spec_Spec_Intercepts.stan")

subsamples <- 25
iterations <- 30
subsampled_species_coef_storage <- array(data = NA, dim = c(iterations, 7, Number_Clusters))

for (i in 1:iterations) {
  set.seed(i)
  
  # Anopheles annularis subsampling
  ann_subsample <- sample(1:Number_Obs_Ann, subsamples, replace = FALSE)
  ann_subsample <- ann_subsample[order(ann_subsample)] 
  Envt_Variables_Ann <- Envt_Variables[species == "Annularis", ]
  Envt_Variables_Ann <- Envt_Variables_Ann[ann_subsample, ] 
  Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 1]
  Cluster_Identity_Ann <- Cluster_Identity_Ann[ann_subsample]
  
  # Anopheles culicifacies subsampling  
  cul_subsample <- sample(1:Number_Obs_Cul, subsamples, replace = FALSE)
  cul_subsample <- cul_subsample[order(cul_subsample)] 
  Envt_Variables_Cul <- Envt_Variables[species == "Culicifacies", ]
  Envt_Variables_Cul <- Envt_Variables_Cul[cul_subsample, ] 
  Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 1]
  Cluster_Identity_Cul <- Cluster_Identity_Cul[cul_subsample]
  
  # Anopheles dirus - no subsampling as these have fewer observations than subsamples
  # dir_subsample <- sample(1:Number_Obs_Dir, subsamples replace = FALSE)
  # dir_subsample <- dir_subsample[order(dir_subsample)] 
  Envt_Variables_Dir <- Envt_Variables[species == "Dirus", ]
  Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 1]
  
  # Anopheles fluviatilis subsampling
  fluv_subsample <- sample(1:Number_Obs_Fluv, subsamples, replace = FALSE)
  fluv_subsample <- fluv_subsample[order(fluv_subsample)] 
  Envt_Variables_Fluv <- Envt_Variables[species == "Fluviatilis", ]
  Envt_Variables_Fluv <- Envt_Variables_Fluv[fluv_subsample, ] 
  Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 1]
  Cluster_Identity_Fluv <- Cluster_Identity_Fluv[fluv_subsample]
  
  # Anopheles minimus - No subsampling as these have fewer observations than subsamples
  # min_subsample <- sample(1:Number_Obs_Min, subsamples, replace = FALSE)
  # min_subsample <- min_subsample[order(min_subsample)] 
  Envt_Variables_Min <- Envt_Variables[species == "Minimus", ]
  Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 1]
  
  # Anopheles stephensi subsampling
  ste_subsample <- sample(1:Number_Obs_Ste, subsamples, replace = FALSE)
  ste_subsample <- ste_subsample[order(ste_subsample)] 
  Envt_Variables_Ste <- Envt_Variables[species == "Stephensi", ]
  Envt_Variables_Ste <- Envt_Variables_Ste[ste_subsample, ] 
  Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 1]
  Cluster_Identity_Ste <- Cluster_Identity_Ste[ste_subsample]

  # Anopheles subpictus subsampling
  sub_subsample <- sample(1:Number_Obs_Sub, subsamples, replace = FALSE)
  sub_subsample <- sub_subsample[order(sub_subsample)] 
  Envt_Variables_Sub <- Envt_Variables[species == "Subpictus", ]
  Envt_Variables_Sub <- Envt_Variables_Sub[sub_subsample, ] 
  Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 1]
  Cluster_Identity_Sub <- Cluster_Identity_Sub[sub_subsample]
  
  # Generating STAN data
  subsampled_STAN_Data <- list(Number_Clusters = Number_Clusters, Number_Envt_Variables = Number_Envt_Variables, 
                               Number_Obs_Ann = min(Number_Obs_Ann, subsamples), 
                               Number_Obs_Cul = min(Number_Obs_Cul, subsamples), 
                               Number_Obs_Dir = min(Number_Obs_Dir, subsamples),
                               Number_Obs_Fluv = min(Number_Obs_Fluv, subsamples), 
                               Number_Obs_Min = min(Number_Obs_Min, subsamples), 
                               Number_Obs_Ste = min(Number_Obs_Ste, subsamples), 
                               Number_Obs_Sub = min(Number_Obs_Sub, subsamples),
                               Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                               Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                               Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                               Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub)
            
  subsampled_STAN_fit <- sampling(STAN_GLM, data = subsampled_STAN_Data, iter = 2500, chains = 3, 
                                  refresh = 500, control = list(max_treedepth = 15, adapt_delta = 0.95))
  
  check_hmc_diagnostics(subsampled_STAN_fit)
  subsampled_fit <- rstan::extract(subsampled_STAN_fit)
  
  saveRDS(subsampled_fit, paste0("Outputs/Logistic_Regression_Output/Subsampling/Subsampled_Multinomial_Logistic_Regression_STAN_Output_", i, ".rds"))
  
  subsampled_spec_coefs <- rbind(apply(subsampled_fit$alpha_ann, 2, mean), apply(subsampled_fit$alpha_cul, 2, mean), apply(subsampled_fit$alpha_dir, 2, mean), 
                                 apply(subsampled_fit$alpha_fluv, 2, mean), apply(subsampled_fit$alpha_min, 2, mean), apply(subsampled_fit$alpha_ste, 2, mean), 
                                 apply(subsampled_fit$alpha_sub, 2, mean))
  subsampled_species_coef_storage[i, , ] <- subsampled_spec_coefs
  
  print(i)
  
}

full_fit <- readRDS("Outputs/Logistic_Regression_Output/Multinomial_Logistic_Regression_STAN_Output.rds")
full_spec_coefs <- rbind(apply(full_fit$alpha_ann, 2, mean), apply(full_fit$alpha_cul, 2, mean), apply(full_fit$alpha_dir, 2, mean), 
                         apply(full_fit$alpha_fluv, 2, mean), apply(full_fit$alpha_min, 2, mean), apply(full_fit$alpha_ste, 2, mean), 
                         apply(full_fit$alpha_sub, 2, mean))

apply(subsampled_species_coef_storage , c(2, 3), mean)

plot(as.vector(full_spec_coefs), as.vector(apply(subsampled_species_coef_storage , c(2, 3), mean)),
     xlab = "full", ylab = "subsampled", xlim = c(-1, 1), ylim = c(-1, 1))
lines(c(-1, 0, 1, 2), c(-1, 0, 1, 2))

full_spec_coefs[4, ]
apply(subsampled_species_coef_storage , c(2, 3), mean)[4, ]








