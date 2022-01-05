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
library(UpSetR); library(dendextend)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
source("Functions/Logistic_Regression_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 

if (prior == "informative") {
  fit <- readRDS("Outputs/Logistic_Regression_Output/Informative_Prior/Reduced_Subset_STAN_Output.rds")
} else if (prior == "uninformative") {
  fit <- readRDS("Outputs/Logistic_Regression_Output/Uninformative_Prior/Reduced_Subset_STAN_Output.rds")
}

#######################################################################################################
##                                                                                                   ##
##                    Analysing and Exploring Logistic Regression Results for SPECIES                ## 
##                                                                                                   ##
#######################################################################################################
spec_coefs <- rbind(apply(fit$alpha_ann, 2, mean), apply(fit$alpha_cul, 2, mean), apply(fit$alpha_dir, 2, mean), 
                    apply(fit$alpha_fluv, 2, mean), apply(fit$alpha_min, 2, mean), apply(fit$alpha_ste, 2, mean), apply(fit$alpha_sub, 2, mean))
row.names(spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
colnames(spec_coefs) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

# Figure 5A - Hierachical Clustering of Species Coefficients
par(mar = c(6, 0, 6, 8))
species_coefficient_clustering <- hclust(dist(spec_coefs))
spec_dend <- species_coefficient_clustering %>%
  as.dendrogram() %>%
  color_branches(k = 7, col = c("#898989", "#2EBEEA",  "#A54D2C", "#29B200" , "#E0521A", "#A948EA", "#E8A50B")) %>%
  color_labels(col = c("#898989", "#2EBEEA",  "#A54D2C", "#29B200" , "#E0521A", "#A948EA", "#E8A50B")) %>%
  set("labels_cex", 1.5) %>%
  set("highlight_branches_lwd", 3)
pdf(file = "Figures/Figure_5/Figure_5A_Species_Coefficient_Dendogram.pdf", width = 6.5, height = 5.5)
par(mar = c(6, 0, 6, 8))
plot(spec_dend, las = 1, ylab = "", xlab = "", main = "", horiz = TRUE, lwd = 3, cex = 2)
dev.off()

