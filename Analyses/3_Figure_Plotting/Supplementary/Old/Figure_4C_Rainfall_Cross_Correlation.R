#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
##                                                                                                   ##
#######################################################################################################
library(zoo); library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); 
library(shinystan); library(ggplot2); library(reshape2); library(deSolve); library(parallel); 
library(matlib); library(matlab); library(pracma);  library(invgamma); library(tictoc); library(numbers); 
library(stringr); library(reshape2); library(cowplot); library(tibble); library(dplyr)
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
source("Functions/Logistic_Regression_Functions.R")
source("Functions/CHIRPS_Rainfall_Processing_Functions.R")
mosquito_data_full <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data_full$Keep
species <- mosquito_data_full$Species
mosquito_data <- as.matrix(mosquito_data_full[, 24:35])
prior <- "informative"
set.seed(58) 

if (prior == "informative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
} else if (prior == "uninformative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Clustering.rds")
}

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

fitted_mosquito_data <- fitted_storage
normalised_mosquito_data <- normalised

#######################################################################################################
##                                                                                                   ##
##         Calculating and Plotting the Cross Correlation Between Mosquito Catch and Rainfall        ##
##                                                                                                   ##
#######################################################################################################
months_length <- c(10, 10, 11, 10, 9, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
leap_year_months_length <- c(10, 10, 11, 10, 10, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
monthly_sum_rainfall_storage <- matrix(nrow = length(mosquito_data[, 1]), ncol = 36)
layout.matrix <- matrix(c(1, 1, 2, 2), nrow = 1, ncol = 4, byrow = TRUE)
layout(mat = layout.matrix)
par(oma = c(4, 4, 4, 4), mar = c(2, 2, 2, 2))
ylim <- c(0, 0.1)

for (i in 1:length(mosquito_data_full$Ref_ID)) {
  single_record_dataframe <- mosquito_data_full[i, ]
  rainfall <- generate_rainfall_vector(single_record_dataframe)
  monthly_rainfall_output <- calculate_monthly_rainfall_totals(single_record_dataframe, rainfall)
  monthly_sum_rainfall_storage[i, ] <- monthly_rainfall_output$sum_monthly
}

mean_monthly_sum <- apply(monthly_sum_rainfall_storage[cluster_labels$Cluster == i, ], 2, mean)
lines(timepoints, normalise_total(mean_monthly_sum) * 100, type = "l", col = palette()[1], lwd = 5)


palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C"))
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
timepoints <- seq(0, 35)
for (i in 1:length(fitted_mosquito_data[, 1])) {
  cluster_label <- cluster_labels$Cluster[i]
  if (i == 1) {
    plot(timepoints, monthly_sum_rainfall_storage[i, ], type = "l", ylim = c(0, 400), xaxt = "n",
         col = adjustcolor(palette()[cluster_label], alpha.f = 0.1), xlab = "Month", ylab = "Rainfall (mm)", las = 1)
    axis(1, at = seq(0, 35, 3), labels = months, las = 1)
  } else {
    lines(timepoints, monthly_sum_rainfall_storage[i, ], type = "l",
          col = adjustcolor(palette()[cluster_label], alpha.f = 0.1))
  }
}
mean_monthly_sum_1 <- apply(monthly_sum_rainfall_storage[cluster_labels$Cluster == 1, ], 2, mean)
lines(timepoints, mean_monthly_sum, type = "l", col = palette()[1], lwd = 5)

mean_monthly_sum_2 <- apply(monthly_sum_rainfall_storage[cluster_labels$Cluster == 2, ], 2, mean)
lines(timepoints, mean_monthly_sum_2, type = "l", col = palette()[2], lwd = 5)

mean_monthly_sum_3 <- apply(monthly_sum_rainfall_storage[cluster_labels$Cluster == 3, ], 2, mean)
lines(timepoints, mean_monthly_sum_3, type = "l", col = palette()[3], lwd = 5)

mean_monthly_sum_4 <- apply(monthly_sum_rainfall_storage[cluster_labels$Cluster == 4, ], 2, mean)
lines(timepoints, mean_monthly_sum_4, type = "l", col = palette()[4], lwd = 5)

ccf_storage <- c()
for (i in 1:length(fitted_mosquito_data[, 1])) {
  cross_corr_call <- ccf(fitted_mosquito_data[i, ], monthly_sum_rainfall_storage[i, ], lag.max = 0, plot = FALSE)
  ccf <- cross_corr_call$acf
  ccf_storage[i] <- ccf
}

cluster_ccfs <- tibble(Cluster = cluster_labels$Cluster, CCF = ccf_storage)
cluster_mean_ccfs <- cluster_ccfs %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(mean = mean(CCF), sd = sd(CCF), se = sd(CCF)/sqrt(n()))
palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C", "#E71D36", "#52AD9C", "#7761B5", "#3F220F", "#D6D84F", "#363537"))
barplot(cluster_mean_ccfs$mean, ylim = c(-0.6, 0.6), las = 1, ylab = "Cross-Correlation With Rainfall",
        col = palette()[1:4], names.arg = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), axis.lty = 1)

par(mar = c(1, 0, 3, 0))
layout.matrix <- matrix(c(1, 2, 
                          3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = layout.matrix)
for (i in 1:4) {
  temp <- normalised_mosquito_data[cluster_labels$Cluster == i, ]
  plot(timepoints, apply(temp, 2, mean) * 100, type = "l", yaxt = "n", lwd = 2, ylim = c(0, 11), col = palette()[i], las = 1, xaxt = "n")
  mean_monthly_sum <- apply(monthly_sum_rainfall_storage[cluster_labels$Cluster == i, ], 2, mean)
  lines(timepoints, normalise_total(mean_monthly_sum) * 100, type = "l", col = palette()[1], lwd = 5)
}

# save PDF as 4 wide and 6 high
