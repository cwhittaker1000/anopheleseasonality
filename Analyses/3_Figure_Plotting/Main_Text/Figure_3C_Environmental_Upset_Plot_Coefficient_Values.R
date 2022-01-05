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

top_x <- function(coefficients, number) {
  index <- order(abs(coefficients), decreasing = TRUE)
  top_coefficients <- coefficients[index[1:number]]
  signed_top_coefficients <- c()
  for (i in 1:length(top_coefficients)) {
    if (sign(top_coefficients[i]) == 1) {
      signed_top_coefficients[i] <- paste0("pos_", names(top_coefficients[i]))
    } else {
      signed_top_coefficients[i] <- paste0("neg_", names(top_coefficients[i]))
    }
  }
  return(signed_top_coefficients)
}

one <- top_x(envt_betas[, 1], 15)
two <- top_x(envt_betas[, 2], 15)
three <- top_x(envt_betas[, 3], 15)
four <- top_x(envt_betas[, 4], 15)

listInput <- list(two = two, three = three, four = four, one = one)
intersections <- list(list("four", "three"),
                      list("two", "three"),
                      list("two", "four"),
                      list("one", "two"), 
                      list("one", "three"), 
                      list("one", "four"))

pdf(file = "Figures/Figure_5/Figure_5B_Environmental_Coefficient_Upset_Plot.pdf", width = 6.5, height = 5.5)
upset(fromList(listInput), keep.order = TRUE, intersections = intersections, set_size.show = FALSE,
      order.by = "freq", point.size = 7.5, line.size = 3.5, 
      mainbar.y.label = "Top Ecological Variable Intersections",
      text.scale = 1.6)
dev.off()