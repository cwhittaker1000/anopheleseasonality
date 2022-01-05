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
number_covariates <- dim(fit$envt_betas[, , 1])[2] 
number_clusters <- 4
envt_betas <- matrix(nrow = number_covariates, ncol = number_clusters)
for (i in 1:number_covariates) {
  envt_betas[i, ] <- apply(fit$envt_betas[ , i, ], 2, mean)
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

pdf(file = "Figures/Figure_4/Figure_4B_Environmental_Coefficient_Heatmap.pdf", width = 6.5, height = 8)
envt_melt <- melt(envt_betas)
envt_melt$Var1 <- factor(envt_melt$Var1, levels = rev(row.names(envt_betas)))
base_size <- 9
ggplot(envt_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red") +
  xlab ("") + 
  ylab("") +
  theme_minimal()
dev.off()

