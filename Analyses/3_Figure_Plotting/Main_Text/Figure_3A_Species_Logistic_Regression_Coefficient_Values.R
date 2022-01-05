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

cluster_spec_melt <- melt(spec_coefs)
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = rev(unique(cluster_spec_melt$Var1[order(cluster_spec_melt$Var1)])))
pdf(file = "Figures/Figure_4/Figure_4A_Species_Coefficient_Heatmap.pdf", width = 6.5, height = 5.5)
ggplot(cluster_spec_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red", 
                       limits = c(-1, 1.5), breaks = c(-0.75, 0, 0.75, 1.5)) + 
  xlab ("") + 
  ylab("") +
  theme_minimal() +
  theme(text = element_text(size = 20))
dev.off()
