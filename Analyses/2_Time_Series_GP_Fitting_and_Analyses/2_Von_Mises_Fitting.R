#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
##                                                                                                   ##
#######################################################################################################
library(MESS); library(numbers); library(factoextra); library(rgl); library(tsne); library(zoo); 
library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); library(shinystan); 
library(ggplot2); library(reshape2); library(deSolve); library(parallel); library(matlib); library(matlab); 
library(pracma); library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(dplyr); 
library(VGAM); library(circular); library(DescTools)
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 


#######################################################################################################
##                                                                                                   ##
##              Loading In Fitted Negative Binomial GP Results and Normalising                       ##
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
##              Von Mises Fitting - Fitting 1 and 2 Component Von Mises Distribution                 ##
##                                                                                                   ##
##    Fits both a 1 and 2 component Von Mises distribution to the normalised time series. Then       ##
##    extracts information from either the 1 or 2 component distribution to characterise the time    ##
##    series.                                                                                        ##
##                                                                                                   ##
#######################################################################################################
storage <- matrix(nrow = length(normalised[, 1]), ncol = 6)
colnames(storage) <- c("1_Mean", "1_K", "2_1_Mean", "2_K", "2_2_Mean", "2_W")
peaks <- c()
mean <- c()
weights <- c()
for (i in 1:length(normalised[, 1])) {
  x <- seq(1, 36, 1)
  data <- normalised[i, ]
  data <- data/AUC(2*pi*x/36, data)
  values <- von_mises_fitting(data, TRUE)
  storage[i, ] <- values
  diff <- abs(storage[i, "2_1_Mean"] - storage[i, "2_2_Mean"])
  weight_temp <- storage[i, "2_W"]
  if ((diff > 2.094395 & diff < 4.18879) & (weight_temp > 0.3 & weight_temp < 0.7)) {
    peaks[i] <- 2
  } else {
    peaks[i] <- 1
  }
  if (peaks[i] == 1) {
    mean[i] <- min(storage[i, "1_Mean"], abs(2*pi - storage[i, "1_Mean"]))
    weights[i] <- 1
  } else if (peaks[i] == 2) {
    mean[i] <- -5
    weights[i] <- max(storage[i, "2_W"], 1 - storage[i, "2_W"])
  }
}

Von_Mises_Properties <- cbind(mean, weights, peaks)
colnames(Von_Mises_Properties) <- c("Mean_Unimodal_Only", "Von_Mises_Comp_Weight", "Number_of_Peaks")
if (prior == "informative") {
  saveRDS(Von_Mises_Properties, file = "Outputs/Von_Mises_Characterisation/Informative_Prior/Informative_Prior_Von_Mises_Properties.rds")
} else if (prior == "uninformative") {
  saveRDS(Von_Mises_Properties, file = "Outputs/Von_Mises_Characterisation/Uninformative_Prior/Uninformative_Prior_Von_Mises_Properties.rds")
}

