#######################################################################################################
##                                                                                                   ##
##                                Initial Loading of Libraries & Data                                ##
##                                                                                                   ##
#######################################################################################################
library(MESS); library(numbers); library(factoextra); library(rgl); library(tsne); library(zoo); 
library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); library(shinystan); 
library(ggplot2); library(reshape2); library(deSolve); library(parallel); library(matlib); library(matlab); 
library(pracma); library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(dplyr); 
library(VGAM); library(rgl); library(car); library(R.utils)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- as.factor(mosquito_data$Species)
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 
clusters <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")

#######################################################################################################
##                                                                                                   ##
##                          Plotting Species Composition of Each Cluster                             ## 
##                                                                                                   ##
#######################################################################################################
colours <- c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C")
species_names <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus") 
total_species_counts <- as.vector(table(species))
names(total_species_counts) <- species_names
par(mfrow = c(2, 2), mar = c(4, 2, 4, 2))

cluster_1 <- species[clusters$Cluster == 1]
table_1 <- table(cluster_1)
names(table_1) <- species_names
x <- barplot(table_1/total_species_counts, ylim = c(0, 1), las = 1, col = colours, xaxt = "n", main = "Cluster 1")
text(x = x - 0.5, y = -0.25, names(table_1), xpd = TRUE, srt = 45)
mtext("A", side = 3, adj = 0, line = 1, font = 2, cex = 2)

cluster_2 <- species[clusters$Cluster == 2]
table_2 <- table(cluster_2)
names(table_2) <- species_names
x <- barplot(table_2/total_species_counts, ylim = c(0, 1), las = 1, col = colours, xaxt = "n", main = "Cluster 2")
text(x = x - 0.5, y = -0.25, names(table_2), xpd = TRUE, srt = 45)
mtext("B", side = 3, adj = 0, line = 1, font = 2, cex = 2)

cluster_3 <- species[clusters$Cluster == 3]
table_3 <- table(cluster_3)
names(table_3) <- species_names
x <- barplot(table_3/total_species_counts, ylim = c(0, 1), las = 1, col = colours, xaxt = "n", main = "Cluster 3")
text(x = x - 0.5, y = -0.25, names(table_3), xpd = TRUE, srt = 45)
mtext("C", side = 3, adj = 0, line = 1, font = 2, cex = 2)

cluster_4 <- species[clusters$Cluster == 4]
table_4 <- table(cluster_4)
names(table_4) <- species_names
x <- barplot(table_4/total_species_counts, ylim = c(0, 1), las = 1, col = colours, xaxt = "n", main = "Cluster 4")
text(x = x - 0.5, y = -0.25, names(table_4), xpd = TRUE, srt = 45)
mtext("D", side = 3, adj = 0, line = 1, font = 2, cex = 2)

