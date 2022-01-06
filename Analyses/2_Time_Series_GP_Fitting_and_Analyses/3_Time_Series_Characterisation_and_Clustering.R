#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
##                                                                                                   ##
#######################################################################################################
library(MESS); library(numbers); library(factoextra); library(rgl); library(tsne); library(zoo); 
library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); library(shinystan); 
library(ggplot2); library(reshape2); library(deSolve); library(parallel); library(matlib); library(matlab); 
library(pracma); library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(dplyr); 
library(VGAM); library(rgl); library(car)
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
loc_id <- mosquito_data$Location_ID
ts_id <- mosquito_data$Time_Series_ID
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 

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


#######################################################################################################
##                                                                                                   ##
##              Characterising each of the time series and their temporal properties                 ## 
##                                                                                                   ##
#######################################################################################################
timepoints <- timepoints_storage[1, ]
num_columns <- 7
num_time_series <- length(normalised[, 1])

# Creating and Filling the Matrix for the Outputs of the Mathematical Operations
output_storage_vector <- matrix(nrow = num_time_series, ncol = num_columns) 
output_storage_vector[, 1] <- apply(normalised, 1, entropic_measure, timepoints = timepoints)
output_storage_vector[, 2] <- apply(normalised, 1, calculate_peak_distance_from_jan, timepoints = timepoints)
output_storage_vector[, 3] <- apply(normalised, 1, points_greater_than_mean_multiple, multiple = 1.6)

if (prior == "informative") {
  Periodic_Kernel_Median_All_Time_Series <- readRDS("Outputs/Negative_Binomial_GP_Fitting/Informative_Prior/Informative_Prior_Median_Periods.rds")
} else if (prior == "uninformative") {
  Periodic_Kernel_Median_All_Time_Series <- readRDS("Outputs/Negative_Binomial_GP_Fitting/Uninformative_Prior/Uninformative_Prior_Median_Periods.rds")
}
Periodic_Kernel_Median <- Periodic_Kernel_Median_All_Time_Series
output_storage_vector[, 4] <- Periodic_Kernel_Median

if (prior == "informative") {
  Von_Mises_Clustering_Factors <- readRDS("Outputs/Von_Mises_Characterisation/Informative_Prior/Informative_Prior_Von_Mises_Properties.rds")
} else if (prior == "uninformative") {
  Von_Mises_Clustering_Factors <- readRDS("Outputs/Von_Mises_Characterisation/Uninformative_Prior/Uninformative_Prior_Von_Mises_Properties.rds")
}
Von_Mises_Clustering_Factors <- Von_Mises_Clustering_Factors
output_storage_vector[, 5:7] <- Von_Mises_Clustering_Factors
cor(output_storage_vector)
output_storage_vector_scaled <- scale(output_storage_vector, center = TRUE, scale = TRUE)
if (prior == "informative") {
  saveRDS(output_storage_vector_scaled, file = "Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Time_Series_Properties.rds")
} else if (prior == "uninformative") {
  saveRDS(output_storage_vector_scaled, file = "Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Time_Series_Properties.rds")
}
colnames(output_storage_vector_scaled) <- c("Entropy", "Dist from Jan", "Mean x1.6", "GP Period", "VM Mean", "VM Weight", "VM Peaks")

#######################################################################################################
##                                                                                                   ##
##                                        Running a PCA                                              ## 
##                                                                                                   ##
#######################################################################################################
PCA <- prcomp(output_storage_vector_scaled)
summary <- summary(PCA)
loadings <- PCA$rotation[, 1:7]
PCA_output <- as.matrix(output_storage_vector_scaled) %*% loadings
prop_exp <- PCA$sdev^2/sum(PCA$sdev^2)
plot(seq(1, 7, 1), prop_exp, ylab = "% Variance Explained", xlab = "PCA Component", las = 1)
lines(seq(1, 7, 1), prop_exp)
points(seq(1, 7, 1), prop_exp, pch = 20)

# Plotting PCA Output By Species
par(mfrow = c(1, 1), mar = c(3, 3, 3, 3))
palette(c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C"))
plot(PCA_output[, 2], PCA_output[, 1], pch = 20, col = as.factor(species), xlab = "PCA Comp 2 (20%)", ylab = "PCA Comp 1 (55%)", cex = 2)
legend(-2.5, 4, c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus"), cex = 0.75, col = palette(), pch = 20)

# Examining Results for Different Numbers of Clusters
#     cluster 1 = culicifacies, cluster 2 = bimodal, cluster 3 = fluviatilis, cluster 4 = all year
par(mfrow = c(1, 1), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(4, 4, 4, 4))
num_clusters <- 4
four_cluster_output <- cluster_characterisation(output_storage_vector_scaled, num_clusters, normalised)
four_clusters <- data.frame("Loc_ID" = loc_id, "TS_ID" = ts_id, "Cluster" = four_cluster_output$Cluster_Numbers)
if (prior == "informative") {
  saveRDS(four_clusters, file = "Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
} else if (prior == "uninformative") {
  saveRDS(four_clusters, file = "Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Clustering.rds")
}

# 3D Plotting 
open3d() 
plot3d(PCA_output[, 1:3], col = four_clusters$Cluster, size = 10, xlab = "", ylab = "", zlab = "", box = FALSE)
for (i in 1:4) {
  cov <- cov(PCA_output[four_clusters$Cluster == i, 1:3])
  mean <- apply(PCA_output[four_clusters$Cluster == i, 1:3], 2, mean)
  plot3d(ellipse3d(cov, centre = mean, level = 0.75, alpha = 0.5), col = palette()[i], alpha = 0.2, add = TRUE)
}
