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
row.names(unregularised_spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
cluster_spec_melt <- melt(unregularised_spec_coefs)
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = rev(unique(cluster_spec_melt$Var1[order(cluster_spec_melt$Var1)])))
ggplot(cluster_spec_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red", 
                       limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4)) + 
  xlab ("") + 
  ylab("")

number_covariates <- dim(unregularised_fit$envt_betas[, , 1])[2] 
number_clusters <- 4
envt_betas <- matrix(nrow = number_covariates, ncol = number_clusters)
for (i in 1:number_covariates) {
  envt_betas[i, ] <- apply(unregularised_fit$envt_betas[ , i, ], 2, mean)
}
temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Mean_Temp_Driest_Quarter")
rain_red <- c("Annual_Rain", "Rain_Seasonality", "CHIRPS_Min", "Rain_Coldest_Quarter")
arid_red <- c("Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
hydro_red <- c("Water_Areas_Occurrence", "Water_Areas_Recurrence", "Flow_Accumulation")
landcover_red <- c("City_Accessibility", "Dominant_Landcover")
reduced_subset <- c(temp_red, rain_red, arid_red, hydro_red, landcover_red)
variable_names <- c(reduced_subset[-length(reduced_subset)], "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8", "LC9", "LC10", "LC11")
row.names(envt_betas) <- variable_names
ordered <- c("Annual_Mean_Temperature", "Mean_Temp_Driest_Quarter", "Temperature_Seasonality",
             "Annual_Rain", "CHIRPS_Min", "Rain_Coldest_Quarter", "Rain_Seasonality",
             "Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD",
             "Flow_Accumulation", "Water_Areas_Occurrence", "Water_Areas_Recurrence", 
             "City_Accessibility", "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8", "LC9", "LC10", "LC11")
envt_betas <- envt_betas[ordered, ]
row.names(envt_betas) <- gsub("_", " ", ordered)
colnames(envt_betas) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

envt_melt <- melt(envt_betas)
envt_melt$Var1 <- factor(envt_melt$Var1, levels = rev(row.names(envt_betas)))
base_size <- 9
ggplot(envt_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red") +
  xlab ("") + 
  ylab("") +
  theme_minimal()


unregularised_spec_coefs <- rbind(apply(unregularised_fit$alpha_ann, 2, mean), apply(unregularised_fit$alpha_cul, 2, mean), apply(unregularised_fit$alpha_dir, 2, mean), 
                                  apply(unregularised_fit$alpha_fluv, 2, mean), apply(unregularised_fit$alpha_min, 2, mean), apply(unregularised_fit$alpha_ste, 2, mean), 
                                  apply(unregularised_fit$alpha_sub, 2, mean))
regularised_fit <- readRDS("Outputs/Logistic_Regression_Output/Multinomial_Logistic_Regression_STAN_Output.rds")
regularised_spec_coefs <- rbind(apply(regularised_fit$alpha_ann, 2, mean), apply(regularised_fit$alpha_cul, 2, mean), apply(regularised_fit$alpha_dir, 2, mean), 
                                apply(regularised_fit$alpha_fluv, 2, mean), apply(regularised_fit$alpha_min, 2, mean), apply(regularised_fit$alpha_ste, 2, mean), 
                                apply(regularised_fit$alpha_sub, 2, mean))

row.names(unregularised_spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
row.names(regularised_spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
colnames(unregularised_spec_coefs) <- colnames(regularised_spec_coefs)  <- c("Clus_1", "Clus_2", "Clus_3", "Clus_4")

unreg_clus1_order <- unregularised_spec_coefs[, 1][order(unregularised_spec_coefs[, 1])]
reg_clus1_order <- regularised_spec_coefs[, 1][order(regularised_spec_coefs[, 1])]

unreg_clus2_order <- unregularised_spec_coefs[, 2][order(unregularised_spec_coefs[, 2])]
reg_clus2_order <- regularised_spec_coefs[, 2][order(regularised_spec_coefs[, 2])]

unreg_clus3_order <- unregularised_spec_coefs[, 3][order(unregularised_spec_coefs[, 3])]
reg_clus3_order <- regularised_spec_coefs[, 3][order(regularised_spec_coefs[, 3])]

unreg_clus4_order <- unregularised_spec_coefs[, 4][order(unregularised_spec_coefs[, 4])]
reg_clus4_order <- regularised_spec_coefs[, 4][order(regularised_spec_coefs[, 4])]

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

#######################################################################################################
##                                                                                                   ##
##              Assessing Predictive Performance of Logistic Regression Output                       ## 
##                                                                                                   ##
#######################################################################################################

# Generating probabilites of belonging to each cluster for each species' time series
reg_pred_prob_Ann <- predict_probabilities(Cluster_Identity_Ann, "Cluster_Probabilities_Ann", regularised_fit, Number_Clusters)
reg_pred_prob_Cul <- predict_probabilities(Cluster_Identity_Cul, "Cluster_Probabilities_Cul", regularised_fit, Number_Clusters)
reg_pred_prob_Dir <- predict_probabilities(Cluster_Identity_Dir, "Cluster_Probabilities_Dir", regularised_fit, Number_Clusters)
reg_pred_prob_Fluv <- predict_probabilities(Cluster_Identity_Fluv, "Cluster_Probabilities_Fluv", regularised_fit, Number_Clusters)
reg_pred_prob_Min <- predict_probabilities(Cluster_Identity_Min, "Cluster_Probabilities_Min", regularised_fit, Number_Clusters)
reg_pred_prob_Ste <- predict_probabilities(Cluster_Identity_Ste, "Cluster_Probabilities_Ste", regularised_fit, Number_Clusters)
reg_pred_prob_Sub <- predict_probabilities(Cluster_Identity_Sub, "Cluster_Probabilities_Sub", regularised_fit, Number_Clusters)

unreg_pred_prob_Ann <- predict_probabilities(Cluster_Identity_Ann, "Cluster_Probabilities_Ann", unregularised_fit, Number_Clusters)
unreg_pred_prob_Cul <- predict_probabilities(Cluster_Identity_Cul, "Cluster_Probabilities_Cul", unregularised_fit, Number_Clusters)
unreg_pred_prob_Dir <- predict_probabilities(Cluster_Identity_Dir, "Cluster_Probabilities_Dir", unregularised_fit, Number_Clusters)
unreg_pred_prob_Fluv <- predict_probabilities(Cluster_Identity_Fluv, "Cluster_Probabilities_Fluv", unregularised_fit, Number_Clusters)
unreg_pred_prob_Min <- predict_probabilities(Cluster_Identity_Min, "Cluster_Probabilities_Min", unregularised_fit, Number_Clusters)
unreg_pred_prob_Ste <- predict_probabilities(Cluster_Identity_Ste, "Cluster_Probabilities_Ste", unregularised_fit, Number_Clusters)
unreg_pred_prob_Sub <- predict_probabilities(Cluster_Identity_Sub, "Cluster_Probabilities_Sub", unregularised_fit, Number_Clusters)

# Generating probabilites of belonging to each cluster for each species' time series
reg_pred_clus_Ann <- predict_cluster(reg_pred_prob_Ann)
reg_pred_clus_Cul <- predict_cluster(reg_pred_prob_Cul)
reg_pred_clus_Dir <- predict_cluster(reg_pred_prob_Dir)
reg_pred_clus_Fluv <- predict_cluster(reg_pred_prob_Fluv)
reg_pred_clus_Min <- predict_cluster(reg_pred_prob_Min)
reg_pred_clus_Ste <- predict_cluster(reg_pred_prob_Ste)
reg_pred_clus_Sub <- predict_cluster(reg_pred_prob_Sub)

unreg_pred_clus_Ann <- predict_cluster(unreg_pred_prob_Ann)
unreg_pred_clus_Cul <- predict_cluster(unreg_pred_prob_Cul)
unreg_pred_clus_Dir <- predict_cluster(unreg_pred_prob_Dir)
unreg_pred_clus_Fluv <- predict_cluster(unreg_pred_prob_Fluv)
unreg_pred_clus_Min <- predict_cluster(unreg_pred_prob_Min)
unreg_pred_clus_Ste <- predict_cluster(unreg_pred_prob_Ste)
unreg_pred_clus_Sub <- predict_cluster(unreg_pred_prob_Sub)

# Calculating number correct
reg_Correct_Ann <- calculate_correct_predictions(Cluster_Identity_Ann, reg_pred_clus_Ann)
reg_Correct_Cul <- calculate_correct_predictions(Cluster_Identity_Cul, reg_pred_clus_Cul)
reg_Correct_Dir <- calculate_correct_predictions(Cluster_Identity_Dir, reg_pred_clus_Dir)
reg_Correct_Fluv <- calculate_correct_predictions(Cluster_Identity_Fluv, reg_pred_clus_Fluv)
reg_Correct_Min <- calculate_correct_predictions(Cluster_Identity_Min, reg_pred_clus_Min)
reg_Correct_Ste <- calculate_correct_predictions(Cluster_Identity_Ste, reg_pred_clus_Ste)
reg_Correct_Sub <- calculate_correct_predictions(Cluster_Identity_Sub, reg_pred_clus_Sub)

unreg_Correct_Ann <- calculate_correct_predictions(Cluster_Identity_Ann, unreg_pred_clus_Ann)
unreg_Correct_Cul <- calculate_correct_predictions(Cluster_Identity_Cul, unreg_pred_clus_Cul)
unreg_Correct_Dir <- calculate_correct_predictions(Cluster_Identity_Dir, unreg_pred_clus_Dir)
unreg_Correct_Fluv <- calculate_correct_predictions(Cluster_Identity_Fluv, unreg_pred_clus_Fluv)
unreg_Correct_Min <- calculate_correct_predictions(Cluster_Identity_Min, unreg_pred_clus_Min)
unreg_Correct_Ste <- calculate_correct_predictions(Cluster_Identity_Ste, unreg_pred_clus_Ste)
unreg_Correct_Sub <- calculate_correct_predictions(Cluster_Identity_Sub, unreg_pred_clus_Sub)

reg_total_correct <- reg_Correct_Ann + reg_Correct_Cul + reg_Correct_Dir + reg_Correct_Fluv + reg_Correct_Min + reg_Correct_Ste + reg_Correct_Sub
unreg_total_correct <- unreg_Correct_Ann + unreg_Correct_Cul + unreg_Correct_Dir + unreg_Correct_Fluv + unreg_Correct_Min + unreg_Correct_Ste + unreg_Correct_Sub
