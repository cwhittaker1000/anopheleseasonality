#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
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
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
loc_id <- mosquito_data$Location_ID
ts_id <- mosquito_data$Time_Series_ID
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 

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

#######################################################################################################
##                                                                                                   ##
##                        Calculating Proportion of Density in X Months                              ##
##                                                                                                   ##
#######################################################################################################
calc_seas <- function(x, months) {
  seas <- c()
  forw <- (months * 3) - 1
  for (i in 1:length(x)) {
    index <- i:(i+forw)
    index[index>36] <- index[index>36] - 36
    temp <- sum(x[index])
    seas <- c(seas, temp)
  }
  return(max(seas))
}

ann <- normalised[species == "Annularis", ]
cul <- normalised[species == "Culicifacies", ]
dir <- normalised[species == "Dirus", ]
fluv <- normalised[species == "Fluviatilis", ]
min <- normalised[species == "Minimus", ]
ste <- normalised[species == "Stephensi", ]
sub <- normalised[species == "Subpictus", ]

ann_mean <- c()
cul_mean <- c()
dir_mean <- c()
fluv_mean <- c()
min_mean <- c()
ste_mean <- c()
sub_mean <- c()

ann_se <- c()
cul_se <- c()
dir_se <- c()
fluv_se <- c()
min_se <- c()
ste_se <- c()
sub_se <- c()

ann_num <- nrow(ann)
cul_num <- nrow(cul)
dir_num <- nrow(dir)
fluv_num <- nrow(fluv)
min_num <- nrow(min)
ste_num <- nrow(ste)
sub_num <- nrow(sub)

for (i in 1:12) {
  ann_mean <- c(ann_mean, mean(apply(ann, 1, calc_seas, i)))
  cul_mean <- c(cul_mean, mean(apply(cul, 1, calc_seas, i)))
  dir_mean <- c(dir_mean, mean(apply(dir, 1, calc_seas, i)))
  fluv_mean <- c(fluv_mean, mean(apply(fluv, 1, calc_seas, i)))
  min_mean <- c(min_mean, mean(apply(min, 1, calc_seas, i)))
  ste_mean <- c(ste_mean, mean(apply(ste, 1, calc_seas, i)))
  sub_mean <- c(sub_mean, mean(apply(sub, 1, calc_seas, i)))
  
  ann_se <- c(ann_se, sd(apply(ann, 1, calc_seas, i))/sqrt(ann_num))
  cul_se <- c(cul_se, sd(apply(cul, 1, calc_seas, i))/sqrt(cul_num))
  dir_se <- c(dir_se, sd(apply(dir, 1, calc_seas, i))/sqrt(dir_num))
  fluv_se <- c(fluv_se, sd(apply(fluv, 1, calc_seas, i))/sqrt(dir_num))
  min_se <- c(min_se, sd(apply(min, 1, calc_seas, i))/sqrt(min_num))
  ste_se <- c(ste_se, sd(apply(ste, 1, calc_seas, i))/sqrt(ste_num))
  sub_se <- c(sub_se, sd(apply(sub, 1, calc_seas, i))/sqrt(sub_num))
}

palette(c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C"))
plot(0:12, c(0, ann_mean), ylim = c(0, 1), type = "l", col = palette()[1], las = 1, xlab = "Number of Months",
     ylab = "Percentage of Total Annual Catch in X Months", lwd = 2)
lines(0:12, c(0, cul_mean), ylim = c(0, 1), type = "l", col = palette()[2], lwd = 2)
lines(0:12, c(0, dir_mean), ylim = c(0, 1), type = "l", col = palette()[3], lwd = 2)
lines(0:12, c(0, fluv_mean), ylim = c(0, 1), type = "l", col = palette()[4], lwd = 2)
lines(0:12, c(0, min_mean), ylim = c(0, 1), type = "l", col = palette()[5], lwd = 2)
lines(0:12, c(0, ste_mean), ylim = c(0, 1), type = "l", col = palette()[6], lwd = 2)
lines(0:12, c(0, sub_mean), ylim = c(0, 1), type = "l", col = palette()[7], lwd = 2)
lines(seq(from = 0, to = 12, length.out = 13), seq(from = 0, to = 1, length.out = 13), 
      col = "black", lty = 3)

layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE)
layout(mat = layout.matrix)

palette(c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C"))
am <- 0.675
plot(0:12, c(0, ann_mean), ylim = c(0, 1), type = "l", col = palette()[1], las = 1, 
     xlab = "Number of Months Considered", ylab = "",
     lwd = 2, yaxt = "n")
axis(side = 4, at = seq(0, 1, 0.2), las = 1)
polygon(c(0:12, rev(0:12)),
        c(0, ann_mean + am * ann_se, rev(ann_mean - am * ann_se), 0),
        col = adjustcolor(palette()[1], alpha.f = 0.1), border = NA)
lines(0:12, c(0, cul_mean), ylim = c(0, 1), type = "l", col = palette()[2], lwd = 2)
polygon(c(0:12, rev(0:12)),
        c(0, cul_mean + am * cul_se, rev(cul_mean - am * cul_se), 0),
        col = adjustcolor(palette()[2], alpha.f = 0.1), border = NA)
lines(0:12, c(0, dir_mean), ylim = c(0, 1), type = "l", col = palette()[3], lwd = 2)
polygon(c(0:12, rev(0:12)),
        c(0, dir_mean + am * dir_se, rev(dir_mean - am * dir_se), 0),
        col = adjustcolor(palette()[3], alpha.f = 0.1), border = NA)
lines(0:12, c(0, fluv_mean), ylim = c(0, 1), type = "l", col = palette()[4], lwd = 2)
polygon(c(0:12, rev(0:12)),
        c(0, fluv_mean + am * fluv_se, rev(fluv_mean - am * fluv_se), 0),
        col = adjustcolor(palette()[4], alpha.f = 0.1), border = NA)
lines(0:12, c(0, min_mean), ylim = c(0, 1), type = "l", col = palette()[5], lwd = 2)
polygon(c(0:12, rev(0:12)),
        c(0, min_mean + am * min_se, rev(min_mean - am * min_se), 0),
        col = adjustcolor(palette()[5], alpha.f = 0.1), border = NA)
lines(0:12, c(0, ste_mean), ylim = c(0, 1), type = "l", col = palette()[6], lwd = 2)
polygon(c(0:12, rev(0:12)),
        c(0, ste_mean + am * ste_se, rev(ste_mean - am * ste_se), 0),
        col = adjustcolor(palette()[6], alpha.f = 0.1), border = NA)
lines(0:12, c(0, sub_mean), ylim = c(0, 1), type = "l", col = palette()[7], lwd = 2)
polygon(c(0:12, rev(0:12)),
        c(0, sub_mean + am * sub_se, rev(sub_mean - am * sub_se), 0),
        col = adjustcolor(palette()[7], alpha.f = 0.1), border = NA)
lines(seq(from = 0, to = 12, length.out = 13), seq(from = 0, to = 1, length.out = 13), 
      col = "black", lty = 3)
lines(rep(4, 13), seq(from = 0, to = 1, length.out = 13), 
      col = "black", lty = 1, lwd = 2)

ann_ind <- apply(ann, 1, calc_seas, 4)
cul_ind <- apply(cul, 1, calc_seas, 4)
dir_ind <- apply(dir, 1, calc_seas, 4)
fluv_ind <- apply(fluv, 1, calc_seas, 4)
min_ind <- apply(min, 1, calc_seas, 4)
ste_ind <- apply(ste, 1, calc_seas, 4)
sub_ind <- apply(sub, 1, calc_seas, 4)

df_seas <- data.frame(species = c(rep("annularis", length(ann_ind)),
                                  rep("culicifacies", length(cul_ind)),
                                  rep("dirus", length(dir_ind)),
                                  rep("fluviatilis", length(fluv_ind)),
                                  rep("minimus", length(min_ind)),
                                  rep("stephensi", length(ste_ind)),
                                  rep("subpictus", length(sub_ind))),
                      seas = c(ann_ind, cul_ind, dir_ind, fluv_ind, min_ind, ste_ind, sub_ind))
df_seas$species <- factor(df_seas$species, levels=c("subpictus", "stephensi", "minimus", "fluviatilis", "dirus", "culicifacies", "annularis"))

par(mfrow = c(1, 1))
set.seed(11)
boxplot(seas ~ species, data = df_seas, border = palette()[7:1], col = NA, outline = FALSE,
        xlab = "", ylab = "% of Total Annual Catch in 4 Months", lty = 1, horizontal = TRUE, las = 1)
stripchart(seas ~ species, data = df_seas, method = "jitter", jitter = 0.25,
           pch = 20, cex = 1.5, col = palette()[7:1], vertical = FALSE, add = TRUE)
df_seas$species <- factor(df_seas$species, levels=c("subpictus", "stephensi", "minimus", "fluviatilis", "dirus", "culicifacies", "annularis"))

par(mfrow = c(1, 1))
df_seas$species <- factor(df_seas$species, levels=rev(c("subpictus", "stephensi", "minimus", "fluviatilis", "dirus", "culicifacies", "annularis")))
set.seed(11)
boxplot(seas ~ species, data = df_seas, border = palette()[1:71], col = NA, outline = FALSE,
        xlab = "", ylab = "% of Total Annual Catch in 4 Months", lty = 1, horizontal = FALSE, las = 1)
stripchart(seas ~ species, data = df_seas, method = "jitter", jitter = 0.25,
           pch = 20, cex = 1.5, col = palette()[1:7], vertical = TRUE, add = TRUE)

t.test(ann_ind, cul_ind)
t.test(ann_ind, dir_ind)
t.test(ann_ind, fluv_ind)
t.test(ann_ind, min_ind)
t.test(ann_ind, ste_ind)
t.test(ann_ind, sub_ind)

t.test(cul_ind, dir_ind)
t.test(cul_ind, fluv_ind)
t.test(cul_ind, min_ind)
t.test(cul_ind, ste_ind)
t.test(cul_ind, sub_ind)

t.test(fluv_ind, dir_ind)
t.test(fluv_ind, min_ind)
t.test(fluv_ind, ste_ind)
t.test(fluv_ind, sub_ind)

t.test(dir_ind, min_ind)
t.test(dir_ind, ste_ind)
t.test(dir_ind, sub_ind)

t.test(min_ind, ste_ind)
t.test(min_ind, sub_ind)

t.test(ste_ind, sub_ind)

######################################################################################################
##                                                                                                   ##
##         Calculating and Plotting the Cross Correlation Between Mosquito Catch and Rainfall        ##
##                                                                                                   ##
#######################################################################################################
# mosquito_data_full <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
# fitted_mosquito_data <- fitted_storage
# normalised_mosquito_data <- normalised
# 
# months_length <- c(10, 10, 11, 10, 9, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
# leap_year_months_length <- c(10, 10, 11, 10, 10, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
# monthly_sum_rainfall_storage <- matrix(nrow = length(mosquito_data[, 1]), ncol = 36)
# layout.matrix <- matrix(c(1, 1, 1, 1,
#                           1, 1, 1, 1,
#                           1, 1, 1, 1, 
#                           2, 2, 2, 2, 
#                           2, 2, 2, 2, 
#                           2, 2, 2, 2, 
#                           3, 4, 5, 6,
#                           3, 4, 5, 6), nrow = 8, ncol = 4, byrow = TRUE)
# layout(mat = layout.matrix)
# par(oma = c(4, 4, 4, 4), mar = c(2, 2, 2, 2))
# ylim <- c(0, 0.1)
# 
# for (i in 1:length(mosquito_data_full$Ref_ID)) {
#   single_record_dataframe <- mosquito_data_full[i, ]
#   rainfall <- generate_rainfall_vector(single_record_dataframe)
#   monthly_rainfall_output <- calculate_monthly_rainfall_totals(single_record_dataframe, rainfall)
#   monthly_sum_rainfall_storage[i, ] <- monthly_rainfall_output$sum_monthly
# }
# 
# 
# # By time-series
# seasonality_mosquito <- apply(normalised_mosquito_data, 1, calc_seas, 3)
# 
# normalised_rainfall <- t(apply(monthly_sum_rainfall_storage, 1, normalise_total))
# seasonality_rainfall <- apply(normalised_rainfall, 1, calc_seas, 3)
# 
# plot(seasonality_rainfall, seasonality_mosquito, ylim = c(0, 1))
