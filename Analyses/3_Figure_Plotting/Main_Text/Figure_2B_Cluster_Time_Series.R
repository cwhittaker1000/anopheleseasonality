#######################################################################################################
##                                                                                                   ##
##                                Initial Loading of Libraries & Data                                ##
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
source("Functions/CHIRPS_Rainfall_Processing_Functions.R")
mosquito_data_full <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data_full$Keep
species <- mosquito_data_full$Species
mosquito_data <- as.matrix(mosquito_data_full[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
set.seed(58) 
clusters <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")

#######################################################################################################
##                                                                                                   ##
##              Loading In Fitted Negative Binomial GP Results, Normalising & Removing               ##
##                                      Low Count Time Series                                        ## 
##                                                                                                   ##
#######################################################################################################
fitted_storage <- matrix(nrow = 272, ncol = 36)
timepoints_storage <- matrix(nrow = 272, ncol = 36)
for (i in 1:272) { 
  temp <- mean_realisation_extract(i, mosquito_data, "informative", FALSE)
  fitted_storage[i, 1:36] <- temp$mean
  timepoints_storage[i, 1:36] <- temp$timepoints
}
normalised <- matrix(nrow = 272, ncol = 36)
for (i in 1:272) {
  temp <-  normalise_total(fitted_storage[i, ])
  normalised[i, 1:length(temp)] <- temp
}

time_series <- normalised

#######################################################################################################
##                                                                                                   ##
##                               Loading In CHIRPS Rainfall Data                                     ##
##                                                                                                   ##
#######################################################################################################
months_length <- c(10, 10, 11, 10, 9, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
leap_year_months_length <- c(10, 10, 11, 10, 10, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
monthly_sum_rainfall_storage <- matrix(nrow = length(mosquito_data[, 1]), ncol = 36)

for (i in 1:length(mosquito_data_full$Ref_ID)) {
  single_record_dataframe <- mosquito_data_full[i, ]
  rainfall <- generate_rainfall_vector(single_record_dataframe)
  monthly_rainfall_output <- calculate_monthly_rainfall_totals(single_record_dataframe, rainfall)
  monthly_sum_rainfall_storage[i, ] <- monthly_rainfall_output$sum_monthly
}

#######################################################################################################
##                                                                                                   ##
##                          Plotting the Time Series By Cluster Label                                ## 
##                                                                                                   ##
#######################################################################################################
palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C", "#E71D36", "#52AD9C", "#7761B5", "#3F220F", "#D6D84F", "#363537"))
clusters <- clusters$Cluster
num_clusters <- length(unique(clusters))
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
timepoints <- seq(0, 11.6666667, length.out = 36)

pdf("Figures/Figure_2/Figure_2B_Cluster_Temporal_Patterns.pdf", height = 6, width = 6.5)
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 1, 1, 5))
max <- 10 
lty <- 2
for (i in 1:num_clusters) {
  
  if (i == 1) {
    cluster <- time_series[clusters == i, ]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    mean_monthly_sum <- apply(monthly_sum_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints, normalise_total(mean_monthly_sum) * 100, type = "l", col = "black", lwd = 1, lty = lty)
  } else if (i == 2) {
    cluster <- time_series[clusters == i, ]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    axis(4, at = seq(0, 10, 2), las = 2)
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
      
    }
    mean_monthly_sum <- apply(monthly_sum_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints, normalise_total(mean_monthly_sum) * 100, type = "l", col = "black", lwd = 1, lty = lty)
  } else if (i == 3) {
    cluster <- time_series[clusters == i, ]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    axis(1, at = seq(0, 11, 1), labels = months, las = 2)
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    mean_monthly_sum <- apply(monthly_sum_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints, normalise_total(mean_monthly_sum) * 100, type = "l", col = "black", lwd = 1, lty = lty)
  } else {
    cluster <- time_series[clusters == i, ]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    axis(4, at = seq(0, 10, 2), las = 2)
    axis(1, at = seq(0, 11, 1), labels = months, las = 2)
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    mean_monthly_sum <- apply(monthly_sum_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints, normalise_total(mean_monthly_sum) * 100, type = "l", col = "black", lwd = 1, lty = lty)
  }
  mtext("Normalised Catch (% of Annual Total)", side = 4, outer = TRUE, cex = 1, font = 2, line = 3, col = "grey20")
}
dev.off()




