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
library(UpSetR); library(dendextend); library(here); library(tidyverse)

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
##                     Running Logistic Regression Models Implemented Using STAN                     ## 
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

Envt_Variables_Ann <- Envt_Variables[species == "Annularis", ]
Envt_Variables_Cul <- Envt_Variables[species == "Culicifacies", ]
Envt_Variables_Dir <- Envt_Variables[species == "Dirus", ]
Envt_Variables_Fluv <- Envt_Variables[species == "Fluviatilis", ]
Envt_Variables_Min <- Envt_Variables[species == "Minimus", ]
Envt_Variables_Ste <- Envt_Variables[species == "Stephensi", ]
Envt_Variables_Sub <- Envt_Variables[species == "Subpictus", ]

Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 1]
Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 1]
Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 1]
Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 1]
Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 1]
Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 1]
Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 1]

STAN_Data <- list(Number_Clusters = Number_Clusters, Number_Envt_Variables = Number_Envt_Variables, 
                  Number_Obs_Ann = Number_Obs_Ann, Number_Obs_Cul = Number_Obs_Cul, Number_Obs_Dir = Number_Obs_Dir,
                  Number_Obs_Fluv = Number_Obs_Fluv, Number_Obs_Min = Number_Obs_Min, Number_Obs_Ste = Number_Obs_Ste, Number_Obs_Sub = Number_Obs_Sub,
                  Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                  Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                  Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                  Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub,
                  tau = 10)

unregularised_STAN_GLM <- stan_model("Model_Files/Multinomial_Spec_Spec_Intercepts_No_Regularisation.stan")
options(mc.cores = 3)
unregularised_STAN_fit <- sampling(unregularised_STAN_GLM, data = STAN_Data, iter = 3000, chains = 3, refresh = 500, control = list(max_treedepth = 15, adapt_delta = 0.95))

check_hmc_diagnostics(unregularised_STAN_fit)
unregularised_fit <- rstan::extract(unregularised_STAN_fit)

saveRDS(unregularised_fit, "Outputs/Logistic_Regression_Output/Unregularised_Multinomial_Logistic_Regression_STAN_Output.rds")

unregularised_spec_coefs <- rbind(apply(unregularised_fit$alpha_ann, 2, mean), apply(unregularised_fit$alpha_cul, 2, mean), apply(unregularised_fit$alpha_dir, 2, mean), 
                                  apply(unregularised_fit$alpha_fluv, 2, mean), apply(unregularised_fit$alpha_min, 2, mean), apply(unregularised_fit$alpha_ste, 2, mean), 
                                  apply(unregularised_fit$alpha_sub, 2, mean))

regularised_fit <- readRDS("Outputs/Logistic_Regression_Output/Multinomial_Logistic_Regression_STAN_Output.rds")
regularised_spec_coefs <- rbind(apply(regularised_fit$alpha_ann, 2, mean), apply(regularised_fit$alpha_cul, 2, mean), apply(regularised_fit$alpha_dir, 2, mean), 
                                apply(regularised_fit$alpha_fluv, 2, mean), apply(regularised_fit$alpha_min, 2, mean), apply(regularised_fit$alpha_ste, 2, mean), 
                                apply(regularised_fit$alpha_sub, 2, mean))

ann_df_unreg <- data.frame(type = "unregularised", unregularised_fit$alpha_ann)
ann_df_reg <- data.frame(type = "regularised", regularised_fit$alpha_ann)
ann_df <- rbind(ann_df_unreg, ann_df_reg) %>%
  pivot_longer(cols = starts_with("X"), names_to = "coef", values_to = "value")


cul_df_unreg <- data.frame(type = "unregularised", unregularised_fit$alpha_cul)
cul_df_reg <- data.frame(type = "regularised", regularised_fit$alpha_cul)
cul_df <- rbind(cul_df_unreg, cul_df_reg) %>%
  pivot_longer(cols = starts_with("X"), names_to = "coef", values_to = "value")

ggplot(cul_df[cul_df$type == "regularised", ], aes(x = value, fill = type)) +
  geom_density() +
  facet_wrap(~coef)

unregularised_spec_coefs/regularised_spec_coefs

apply(abs(unregularised_spec_coefs), 1, sum)
apply(abs(regularised_spec_coefs), 1, sum)
