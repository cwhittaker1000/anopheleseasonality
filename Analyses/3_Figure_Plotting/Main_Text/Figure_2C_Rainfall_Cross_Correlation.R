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

for (i in 1:length(mosquito_data_full$Ref_ID)) {
  single_record_dataframe <- mosquito_data_full[i, ]
  rainfall <- generate_rainfall_vector(single_record_dataframe)
  monthly_rainfall_output <- calculate_monthly_rainfall_totals(single_record_dataframe, rainfall)
  monthly_sum_rainfall_storage[i, ] <- monthly_rainfall_output$sum_monthly
}

pdf("Figures/Figure_2/Figure_2C_Rainfall_Cross_Correlations.pdf", height = 5.5, width = 5.5)

palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C"))
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

par(mfrow = c(1, 1))
cluster_ccfs$Cluster <- factor(cluster_ccfs$Cluster, levels=c(4, 3, 2, 1))

boxplot(CCF ~ Cluster, data = cluster_ccfs, border = palette()[4:1], col = NA, horizontal = TRUE,
        lex.order = FALSE, names = c("Cluster 4", "Cluster 3", "Cluster 2", "Cluster 1"), outline = FALSE, las = 1,
        ylim = c(-1, 1), ylab = "", xlab = "")
stripchart(CCF ~ Cluster, data = cluster_ccfs, method = "jitter", jitter = 0.25,
           pch = 20, cex = 1.5, col = palette()[4:1], vertical = FALSE, add = TRUE)
dev.off()