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
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Informative_Prior/Informative_Prior_Clustering.rds")
} else if (prior == "uninformative") {
  cluster_labels <- readRDS("Outputs/Characterisation_and_Clustering/Uninformative_Prior/Uninformative_Prior_Clustering.rds")
}
covariates <- readRDS("Datasets/Extracted_Covariates_for_Modelling/Extracted_Covariates.rds")


#######################################################################################################
##                                                                                                   ##
##        Loading In Environmental Covariates and Cluster Assignments for Each Time Series           ## 
##                                                                                                   ##
#######################################################################################################

# Reduced Subset (4 + 7 + 4 + 6 + 2 + 11 landcover variables) - 34 variables
temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Day_LST_SD", "Mean_Temp_Driest_Quartest")
rain_red <- c("Annual_Rain", "Rain_Seasonality", "WC_A3", "WC_P0", "WC_P3", "CHIRPS_Min", "Rain_Coldest_Quarter")
arid_red <- c("PET_Yearly_Average", "Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
hydro_red <- c("WWF_Distance_to_Water", "Water_Areas_Recurrence", "Water_Areas_Seasonality", "Irrigated_Areas", "Flow_Accumulation", "Elevation")
landcover_red <- c("City_Accessibility", "EVI_Mean", "Dominant_Landcover")
reduced_subset <- c(temp_red, rain_red, arid_red, hydro_red, landcover_red)

# Even More Reduced Subset (3 + 4 + 3 + 3 + 1 + 11 landcover variables) - 25 variables
temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Mean_Temp_Driest_Quartest")
rain_red <- c("Annual_Rain", "Rain_Seasonality", "CHIRPS_Min", "Rain_Coldest_Quarter")
arid_red <- c("Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
hydro_red <- c("Water_Areas_Occurrence", "Water_Areas_Recurrence", "Flow_Accumulation")
landcover_red <- c("City_Accessibility", "Dominant_Landcover")
reduced_subset <- c(temp_red, rain_red, arid_red, hydro_red, landcover_red)

covariates <- covariates[order(covariates[, "location_IDs"]), ]
covariates <- covariates[, -1]
covariates <- as.data.frame(covariates)
covariates <- covariates[, reduced_subset]
covariates$Dominant_Landcover <- as.factor(covariates$Dominant_Landcover)
covariates <- as.matrix(one_hot(as.data.table(covariates), cols = "Dominant_Landcover"))
number_covariates <- dim(covariates)[2] 

correlation <- melt(cor(covariates))
ggplot(correlation, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

number_time_series <- length(cluster_labels$Cluster)
temp <- cluster_labels$Cluster
temp <- as.data.frame(cbind(Overall_Cluster = temp, Cluster = temp))
temp$Cluster <- as.factor(temp$Cluster)
Cluster_Labels <- as.data.frame(one_hot(as.data.table(temp), cols = "Cluster"))

# Creating the Datasets
Envt_Variables <- matrix(nrow = number_time_series, ncol = number_covariates)
for (i in 1:number_time_series) {
  index <- cluster_labels$Loc_ID[i]
  Envt_Variables[i, ] <- unname(as.vector(covariates[index, ]))
  print(i)
}
variable_names <- c(reduced_subset[-length(reduced_subset)], "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8", "LC9", "LC10", "LC11")
colnames(Envt_Variables) <- variable_names


#######################################################################################################
##                                                                                                   ##
##                     Running Logistic Regression Models Implemented Using STAN                     ## 
##                                                                                                   ##
#######################################################################################################
Number_Envt_Variables <- number_covariates

Number_Obs_Ann <- sum(species == "Annularis")
Number_Obs_Cul <- sum(species == "Culicifacies")
Number_Obs_Dir <- sum(species == "Dirus")
Number_Obs_Fluv <- sum(species == "Fluviatilis")
Number_Obs_Min <- sum(species == "Minimus")
Number_Obs_Ste <- sum(species == "Stephensi")
Number_Obs_Sub <- sum(species == "Subpictus")

Envt_Variables_Ann <- Envt_Variables[species == "Annularis", ]
Envt_Variables_Cul <- Envt_Variables[species == "Culicifacies", ]
Envt_Variables_Dir <- Envt_Variables[species == "Dirus", ]
Envt_Variables_Fluv <- Envt_Variables[species == "Fluviatilis", ]
Envt_Variables_Min <- Envt_Variables[species == "Minimus", ]
Envt_Variables_Ste <- Envt_Variables[species == "Stephensi", ]
Envt_Variables_Sub <- Envt_Variables[species == "Subpictus", ]

STAN_GLM <- stan_model("Model_Files/Logistic_Regression_Spec_Spec_Intercepts.stan")
options(mc.cores = 4)

Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 2]
Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 2]
Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 2]
Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 2]
Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 2]
Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 2]
Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 2]

STAN_Data <- list(Number_Envt_Variables = Number_Envt_Variables, 
                  Number_Obs_Ann = Number_Obs_Ann, Number_Obs_Cul = Number_Obs_Cul, Number_Obs_Dir = Number_Obs_Dir,
                  Number_Obs_Fluv = Number_Obs_Fluv, Number_Obs_Min = Number_Obs_Min, Number_Obs_Ste = Number_Obs_Ste, Number_Obs_Sub = Number_Obs_Sub,
                  Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                  Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                  Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                  Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub)

STAN_fit <- sampling(STAN_GLM, data = STAN_Data, iter = 10000, chains = 2, refresh = 5000, control = list(max_treedepth = 15, adapt_delta = 0.95))
fit_1 <- rstan::extract(STAN_fit)

Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 3]
Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 3]
Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 3]
Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 3]
Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 3]
Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 3]
Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 3]

STAN_Data <- list(Number_Envt_Variables = Number_Envt_Variables, 
                  Number_Obs_Ann = Number_Obs_Ann, Number_Obs_Cul = Number_Obs_Cul, Number_Obs_Dir = Number_Obs_Dir,
                  Number_Obs_Fluv = Number_Obs_Fluv, Number_Obs_Min = Number_Obs_Min, Number_Obs_Ste = Number_Obs_Ste, Number_Obs_Sub = Number_Obs_Sub,
                  Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                  Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                  Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                  Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub)

STAN_fit <- sampling(STAN_GLM, data = STAN_Data, iter = 10000, chains = 2, refresh = 5000, control = list(max_treedepth = 15, adapt_delta = 0.95))
fit_2 <- rstan::extract(STAN_fit)


Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 4]
Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 4]
Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 4]
Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 4]
Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 4]
Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 4]
Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 4]

STAN_Data <- list(Number_Envt_Variables = Number_Envt_Variables, 
                  Number_Obs_Ann = Number_Obs_Ann, Number_Obs_Cul = Number_Obs_Cul, Number_Obs_Dir = Number_Obs_Dir,
                  Number_Obs_Fluv = Number_Obs_Fluv, Number_Obs_Min = Number_Obs_Min, Number_Obs_Ste = Number_Obs_Ste, Number_Obs_Sub = Number_Obs_Sub,
                  Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                  Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                  Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                  Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub)

STAN_fit <- sampling(STAN_GLM, data = STAN_Data, iter = 10000, chains = 2, refresh = 5000, control = list(max_treedepth = 15, adapt_delta = 0.95))
fit_3 <- rstan::extract(STAN_fit)


Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 5]
Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 5]
Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 5]
Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 5]
Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 5]
Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 5]
Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 5]

STAN_Data <- list(Number_Envt_Variables = Number_Envt_Variables, 
                  Number_Obs_Ann = Number_Obs_Ann, Number_Obs_Cul = Number_Obs_Cul, Number_Obs_Dir = Number_Obs_Dir,
                  Number_Obs_Fluv = Number_Obs_Fluv, Number_Obs_Min = Number_Obs_Min, Number_Obs_Ste = Number_Obs_Ste, Number_Obs_Sub = Number_Obs_Sub,
                  Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                  Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                  Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                  Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub)

STAN_fit <- sampling(STAN_GLM, data = STAN_Data, iter = 10000, chains = 2, refresh = 5000, control = list(max_treedepth = 15, adapt_delta = 0.95))
fit_4 <- rstan::extract(STAN_fit)


#######################################################################################################
##                                                                                                   ##
##                         Analysing and Exploring Logistic Regression Results                       ## 
##                                                                                                   ##
#######################################################################################################
spec_coefs_1 <- rbind(mean(fit_1$alpha_ann), mean(fit_1$alpha_cul), mean(fit_1$alpha_dir), 
                    mean(fit_1$alpha_fluv), mean(fit_1$alpha_min), mean(fit_1$alpha_ste), mean(fit_1$alpha_sub))
spec_coefs_2 <- rbind(mean(fit_2$alpha_ann), mean(fit_2$alpha_cul), mean(fit_2$alpha_dir), 
                      mean(fit_2$alpha_fluv), mean(fit_2$alpha_min), mean(fit_2$alpha_ste), mean(fit_2$alpha_sub))
spec_coefs_3 <- rbind(mean(fit_3$alpha_ann), mean(fit_3$alpha_cul), mean(fit_3$alpha_dir), 
                      mean(fit_3$alpha_fluv), mean(fit_3$alpha_min), mean(fit_3$alpha_ste), mean(fit_3$alpha_sub))
spec_coefs_4 <- rbind(mean(fit_4$alpha_ann), mean(fit_4$alpha_cul), mean(fit_4$alpha_dir), 
                      mean(fit_4$alpha_fluv), mean(fit_4$alpha_min), mean(fit_4$alpha_ste), mean(fit_4$alpha_sub))
spec_coefs <- cbind(spec_coefs_1, spec_coefs_2, spec_coefs_3, spec_coefs_4)
colnames(spec_coefs) <- c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4")
row.names(spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
cluster_spec_melt <- melt(spec_coefs)
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = rev(unique(cluster_spec_melt$Var1[order(cluster_spec_melt$Var1)])))
ggplot(cluster_spec_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")

spec_melt_cor <- melt(cor(spec_coefs))
base_size <- 9
ggplot(spec_melt_cor, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")

envt_betas <- matrix(nrow = number_covariates, ncol = 4)
for (i in 1:number_covariates) {
  for (j in 1:4) {
    if (j == 1) {
      envt_betas[i, j] <- mean(fit_1$envt_betas[ , i])
    } else if (j == 2) {
      envt_betas[i, j] <- mean(fit_2$envt_betas[ , i])
    } else if (j == 3) {
      envt_betas[i, j] <- mean(fit_3$envt_betas[ , i])
    } else if (j == 4) {
      envt_betas[i, j] <- mean(fit_4$envt_betas[ , i])
    }
  }

}
row.names(envt_betas) <- variable_names
colnames(envt_betas) <- c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4")
envt_melt <- melt(envt_betas)
envt_melt$Var1 <- factor(envt_melt$Var1, levels = rev(unique(envt_melt$Var1[order(envt_melt$Var1)])))
base_size <- 9
ggplot(envt_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")

cluster_1_betas <- envt_betas[, 1]
index_for_ordering <- order(cluster_1_betas)
envt_melt_test <- melt(envt_betas)
envt_melt_test$Var1 <- factor(envt_melt_test$Var1, levels = envt_melt_test$Var1[index_for_ordering])
base_size <- 9
ggplot(envt_melt_test, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")

envt_melt_cor <- melt(cor(envt_betas))
envt_melt_cor$Var1 <- factor(envt_melt_cor$Var1, levels = rev(unique(envt_melt_cor$Var1[order(envt_melt_cor$Var1)])))
base_size <- 9
ggplot(envt_melt_cor, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")


#######################################################################################################
##                                                                                                   ##
##              Analysing and Exploring Logistic Regression Results for ENVIRONMENTAL                ## 
##                                                                                                   ##
#######################################################################################################
one_pos <- paste("pos", variable_names[order(envt_betas[, 1], decreasing = TRUE)[1:10]])
one_neg <- paste("neg", variable_names[order(envt_betas[, 1], decreasing = FALSE)[1:10]])
one <- c(one_pos, one_neg)

two_pos <- paste("pos", variable_names[order(envt_betas[, 2], decreasing = TRUE)[1:10]])
two_neg <- paste("neg", variable_names[order(envt_betas[, 2], decreasing = FALSE)[1:10]])
two <- c(two_pos, two_neg)

three_pos <- paste("pos", variable_names[order(envt_betas[, 3], decreasing = TRUE)[1:10]])
three_neg <- paste("neg", variable_names[order(envt_betas[, 3], decreasing = FALSE)[1:10]])
three <- c(three_pos, three_neg)

four_pos <- paste("pos", variable_names[order(envt_betas[, 4], decreasing = TRUE)[1:10]])
four_neg <- paste("neg", variable_names[order(envt_betas[, 4], decreasing = FALSE)[1:10]])
four <- c(four_pos, four_neg)

listInput <- list(two = two, three = three, four = four, one = one)
intersections <- list(list("four", "three"),
                      list("two", "three"),
                      list("two", "four"),
                      list("one", "two"), 
                      list("one", "three"), 
                      list("one", "four"))
upset(fromList(listInput), keep.order = TRUE, intersections = intersections, set_size.show = FALSE,
      order.by = "freq", point.size = 7.5, line.size = 3.5, 
      mainbar.y.label = "Top Ecological Variable Intersections",
      text.scale = 1.6)
