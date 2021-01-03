#######################################################################################################
##                                                                                                   ##
##                                Initial Loading of Libraries & Data                                ##
##                                                                                                   ##
#######################################################################################################
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
mosquito_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"


#######################################################################################################
##                                                                                                   ##
##                          Loading In Raw Mosquito Data and Plotting Them                           ## 
##                                                                                                   ##
#######################################################################################################
normalised <- matrix(nrow = length(mosquito_data[, 1]), ncol = 12)
for (i in 1:length(mosquito_data[, 1])) {
  temp <-  normalise_total(mosquito_data[i, ])
  normalised[i, 1:length(temp)] <- temp
}

layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4, byrow = TRUE)
layout(mat = layout.matrix)
par(oma = c(3, 6, 2, 3), mar = c(0.5, 0.5, 0.5, 0.5))
ylim <- c(0, 30)
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
colours <- c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C")
timepoints <- seq(1, 12)

annularis <- normalised[species == "Annularis", ]
n_ann <- dim(annularis)[1]
mean_ann <- apply(annularis, 2, mean) * 100
plot(timepoints, annularis[1, ] * 100, type = "l", ylim = ylim, yaxt = "n", las = 1, xaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[1], alpha.f = 0.1))
text(6, 29, paste0("Anopheles annularis"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_ann, ")"), cex = 1.25, font = 2)
axis(2, at = seq(0, 30, 5), las = 2)
for (i in 2:length(annularis[, 1])) {
  lines(timepoints, annularis[i, ] * 100, type = "l", col = adjustcolor(colours[1], alpha.f = 0.1))
}
lines(timepoints, mean_ann, lwd = 3, col = colours[1])

culicifacies <- normalised[species == "Culicifacies", ]
n_cul <- dim(culicifacies)[1]
mean_cul <- apply(culicifacies, 2, mean) * 100
plot(culicifacies[1, ] * 100, type = "l", ylim = ylim, las = 1, xaxt = "n",  yaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[2], alpha.f = 0.1))
text(6, 29, paste0("Anopheles culicifacies"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_cul, ")"), cex = 1.25, font = 2)
for (i in 2:length(culicifacies[, 1])) {
  lines(culicifacies[i, ] * 100, type = "l", col = adjustcolor(colours[2], alpha.f = 0.1))
}
lines(mean_cul, lwd = 3, col = colours[2])

dirus <- normalised[species == "Dirus", ]
n_dir <- dim(dirus)[1]
mean_dirus <- apply(dirus, 2, mean) * 100
plot(dirus[1, ] * 100, type = "l", ylim = ylim, las = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[3], alpha.f = 0.1))
text(6, 29, paste0("Anopheles dirus"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_dir, ")"), cex = 1.25, font = 2)
for (i in 2:length(dirus[, 1])) {
  lines(dirus[i, ] * 100, type = "l", col = adjustcolor(colours[3], alpha.f = 0.1))
}
lines(mean_dirus, lwd = 3, col = colours[3])

fluviatilis <- normalised[species == "Fluviatilis", ]
n_fluv <- dim(fluviatilis)[1]
mean_fluv <- apply(fluviatilis, 2, mean) * 100
plot(fluviatilis[1, ] * 100, type = "l", ylim = ylim, las = 1, xaxt = "n",  yaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[4], alpha.f = 0.1))
text(6, 29, paste0("Anopheles fluviatilis"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_fluv, ")"), cex = 1.25, font = 2)
axis(1, at = seq(1, 12, 1), labels = months, las = 2)
for (i in 2:length(fluviatilis[, 1])) {
  lines(fluviatilis[i, ] * 100, type = "l", col = adjustcolor(colours[4], alpha.f = 0.1))
}
lines(mean_fluv, lwd = 3, col = colours[4])

minimus <- normalised[species == "Minimus", ]
n_min <- dim(minimus)[1]
mean_minimus <- apply(minimus, 2, mean) * 100
plot(minimus[1, ] * 100, type = "l", ylim = ylim, las = 1, yaxt = "n", xaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[5], alpha.f = 0.1))
axis(1, at = seq(1, 12, 1), labels = months, las = 2)
axis(2, at = seq(0, 30, 5), las = 2)
text(6, 29, paste0("Anopheles minimus"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_min, ")"), cex = 1.25, font = 2)
for (i in 2:length(minimus[, 1])) {
  lines(minimus[i, ] * 100, type = "l", col = adjustcolor(colours[5], alpha.f = 0.1))
}
lines(mean_minimus, lwd = 3, col = colours[5])

stephensi <- normalised[species == "Stephensi", ]
n_ste <- dim(stephensi)[1]
mean_steph <- apply(stephensi, 2, mean) * 100
plot(stephensi[1, ] * 100, type = "l", ylim = ylim, las = 1, xaxt = "n",  yaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[6], alpha.f = 0.1))
text(6, 29, paste0("Anopheles stephensi"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_ste, ")"), cex = 1.25, font = 2)
axis(1, at = seq(1, 12, 1), labels = months, las = 2)
for (i in 2:length(stephensi[, 1])) {
  lines(stephensi[i, ] * 100, type = "l", col = adjustcolor(colours[6], alpha.f = 0.1))
}
lines(mean_steph, lwd = 3, col = colours[6])

subpictus <- normalised[species == "Subpictus", ]
n_sub <- dim(subpictus)[1]
mean_sub <- apply(subpictus, 2, mean) * 100
plot(subpictus[1, ] * 100, type = "l", ylim = ylim, las = 1, xaxt = "n",  yaxt = "n", xlab = "", ylab = "", col = adjustcolor(colours[7], alpha.f = 0.1))
axis(1, at = seq(1, 12, 1), labels = months, las = 2)
text(6, 29, paste0("Anopheles subpictus"), cex = 1.25, font = 4)
text(6, 26.5, paste0("(n = ", n_sub, ")"), cex = 1.25, font = 2)
for (i in 2:length(subpictus[, 1])) {
  lines(subpictus[i, ] * 100, type = "l", col = adjustcolor(colours[7], alpha.f = 0.1))
}
lines(mean_sub, lwd = 3, col = colours[7])

mtext("Normalised Catch (% of Annual Total)", side = 2, outer = TRUE, cex = 1.25, font = 2, line = 3, col = "grey20")

