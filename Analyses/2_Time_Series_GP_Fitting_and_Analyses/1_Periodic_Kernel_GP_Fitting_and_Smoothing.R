#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
##                                                                                                   ##
#######################################################################################################
library(zoo); library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); 
library(shinystan); library(ggplot2); library(reshape2); library(deSolve); library(parallel); 
library(matlib); library(matlab); library(pracma); library(rstan); library(ggplot2); library(invgamma); 
library(tictoc)
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"


#######################################################################################################
##                                                                                                   ##
##                            Negative Binomial Gaussian Process Fitting                             ##
##                                                                                                   ##
##    This next section fits a Negative Binomial Gaussian Process with Periodic Kernel to each of    ##
##    the time series, implemented in the probabilistic programming language STAN.                   ##
##                                                                                                   ##
##    Function used for fitting gives you the option of using either Informative or                  ##
##    Uninformative Priors                                                                           ##
##                                                                                                   ##
#######################################################################################################
options(mc.cores = parallel::detectCores() - 4)
GP_model <- stan_model("Model_Files/Neg_Binom_GP_Periodic_Kernel.stan")
par(mfrow = c(1, 1))

number_interpolating_points <- 2
counter <- 1
median_periods <- c()
for (i in counter:272) {
  input_time_series <- mosquito_data[i, ]
  fitting_output <- Periodic_NegBinom_GP_Fitting(input_time_series, number_interpolating_points, i, prior, TRUE)
  fitting_output <- fitting_output$chain_output
  median_periods[counter] <- median(fitting_output[, "period"])
  counter <- counter + 1
  print(i)
}
if (prior == "informative") {
  saveRDS(median_periods, file = "Outputs/Negative_Binomial_GP_Fitting/Informative_Prior/Informative_Prior_Median_Periods.rds")
} else if (prior == "uninformative") {
  saveRDS(median_periods, file = "Outputs/Negative_Binomial_GP_Fitting/Uninformative_Prior/Uninformative_Prior_Median_Periods.rds")
}

# Visualising the Outputs Manually 
for (i in 1:272) {
  mean_realisation_extract(i, mosquito_data, prior, TRUE)
  text(0, 0, paste0(i, " ", i))
  browser()
}


