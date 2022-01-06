#######################################################################################################
##                                                                                                   ##
##                   Loading Required Libraries and Processing Catch Data                            ##
##                                                                                                   ##
#######################################################################################################
library(MESS); library(numbers); library(factoextra); library(rgl); library(tsne); library(zoo); 
library(forecast); library(TSA); library(mgcv); library(GPfit); library(rstan); library(shinystan); 
library(ggplot2); library(reshape2); library(deSolve); library(parallel); library(matlib); library(matlab); 
library(pracma); library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(dplyr); 
library(VGAM); library(rgl); library(car); library(tidyverse)
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
loc_id <- mosquito_data$Location_ID
ts_id <- mosquito_data$Time_Series_ID
keep <- mosquito_data$Keep
species <- mosquito_data$Species

mosquito_catch_data <- as.matrix(mosquito_data[, 24:35])
colnames(mosquito_data) <- seq(1, 12)
prior <- "informative"
set.seed(58) 

if (prior == "informative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
} else if (prior == "uninformative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Clustering.rds")
}

df <- data.frame(Cluster = factor(cluster_labels$Cluster), catch = summed_catch)
ggplot(df, aes(x = log(catch), fill = Cluster)) +
  geom_density(alpha = 0.5) +
  xlab("Log(Total Mosquito Catch Per Study)") +
  ylab("Density") +
  scale_fill_manual(values = c("#F15025", "#7ACC70", "#00A7E1", "#F2328C")) +
  theme_bw()
ggsave("Figures/Supp_Figures/Supp_Figure_6_Catch_Size_Archetype_Membership.pdf",
       width = 6.5, height = 4.5, units = "in")

summed_catch <- apply(mosquito_catch_data, 1, sum)
median(summed_catch[cluster_labels$Cluster == 1])
median(summed_catch[cluster_labels$Cluster == 2])
median(summed_catch[cluster_labels$Cluster == 3])
median(summed_catch[cluster_labels$Cluster == 4])

mood.test(summed_catch[cluster_labels$Cluster == 1], summed_catch[cluster_labels$Cluster == 2])
mood.test(summed_catch[cluster_labels$Cluster == 1], summed_catch[cluster_labels$Cluster == 3])
mood.test(summed_catch[cluster_labels$Cluster == 1], summed_catch[cluster_labels$Cluster == 4])
mood.test(summed_catch[cluster_labels$Cluster == 2], summed_catch[cluster_labels$Cluster == 3])
mood.test(summed_catch[cluster_labels$Cluster == 2], summed_catch[cluster_labels$Cluster == 4])
mood.test(summed_catch[cluster_labels$Cluster == 3], summed_catch[cluster_labels$Cluster == 4])

mean(summed_catch[cluster_labels$Cluster == 1])
mean(summed_catch[cluster_labels$Cluster == 2])
mean(summed_catch[cluster_labels$Cluster == 3])
mean(summed_catch[cluster_labels$Cluster == 4])

sd(summed_catch[cluster_labels$Cluster == 1])
sd(summed_catch[cluster_labels$Cluster == 2])
sd(summed_catch[cluster_labels$Cluster == 3])
sd(summed_catch[cluster_labels$Cluster == 4])

t.test(summed_catch[cluster_labels$Cluster == 1], summed_catch[cluster_labels$Cluster == 2])
t.test(summed_catch[cluster_labels$Cluster == 1], summed_catch[cluster_labels$Cluster == 3])
t.test(summed_catch[cluster_labels$Cluster == 1], summed_catch[cluster_labels$Cluster == 4])
t.test(summed_catch[cluster_labels$Cluster == 2], summed_catch[cluster_labels$Cluster == 3])
t.test(summed_catch[cluster_labels$Cluster == 2], summed_catch[cluster_labels$Cluster == 4])
t.test(summed_catch[cluster_labels$Cluster == 3], summed_catch[cluster_labels$Cluster == 4])
