#######################################################################################################
##                                                                                                   ##
##                                Initial Loading of Libraries & Data                                ##
##                                                                                                   ##
#######################################################################################################
library(MESS); library(numbers); library(factoextra); library(rgl); library(tsne); library(zoo); 
library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); library(shinystan); 
library(ggplot2); library(reshape2); library(deSolve); library(parallel); library(matlib); 
library(matlab); library(pracma); library(rstan); library(ggplot2); library(invgamma); library(tictoc); 
library(dplyr); library(VGAM); library(rgl); library(car)
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
set.seed(58) 
time_series_properties <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Time_Series_Properties.rds")
clusters <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")

#######################################################################################################
##                                                                                                   ##
##                                        Running a PCA                                              ## 
##                                                                                                   ##
#######################################################################################################
PCA <- prcomp(time_series_properties)
summary(PCA)
loadings <- PCA$rotation[, 1:7]
PCA_output <- as.matrix(time_series_properties) %*% loadings
clusters <- clusters$Cluster

# 3D Plotting
palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C", "#E71D36", "#52AD9C", "#7761B5", "#3F220F", "#D6D84F", "#363537"))
open3d() 
plot3d(PCA_output[, 1:3], col = clusters, size = 10, xlab = "", ylab = "", zlab = "", box = FALSE)
for (i in 1:4) {
  cov <- cov(PCA_output[clusters == i, 1:3])
  mean <- apply(PCA_output[clusters == i, 1:3], 2, mean)
  plot3d(ellipse3d(cov, centre = mean, level = 0.75, alpha = 0.5), col = palette()[i], alpha = 0.2, add = TRUE)
}
variance_explained <- summary(PCA)$importance
write.csv(variance_explained, file = "Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Four_Clusters_Reduced_Subset_PCA_Proportion_Variance_Explained.png")
rgl.snapshot(filename = "Figures/Figure_2/Figure_2A_PCA_Output.png")

