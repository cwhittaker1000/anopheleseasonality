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
library(UpSetR); library(dendextend); library(here)

# Sourcing functions
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
source("Functions/Logistic_Regression_Functions.R")

# Loading in and processing mosquito data
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
# Even More Reduced Subset (3 + 4 + 3 + 3 + 1 + 11 landcover variables) - 25 variables
temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Mean_Temp_Driest_Quarter")
rain_red <- c("Annual_Rain", "Rain_Seasonality", "CHIRPS_Min", "Rain_Coldest_Quarter")
arid_red <- c("Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
hydro_red <- c("Water_Areas_Occurrence", "Water_Areas_Recurrence", "Flow_Accumulation")
landcover_red <- c("City_Accessibility", "Dominant_Landcover")
reduced_subset <- c(temp_red, rain_red, arid_red, hydro_red, landcover_red)

colnames(covariates)[10] <- "Mean_Temp_Driest_Quarter"
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
Number_Clusters <- length(unique(cluster_labels$Cluster))
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

Cluster_Identity_Ann <- Cluster_Labels[species == "Annularis", 1]
Cluster_Identity_Cul <- Cluster_Labels[species == "Culicifacies", 1]
Cluster_Identity_Dir <- Cluster_Labels[species == "Dirus", 1]
Cluster_Identity_Fluv <- Cluster_Labels[species == "Fluviatilis", 1]
Cluster_Identity_Min <- Cluster_Labels[species == "Minimus", 1]
Cluster_Identity_Ste <- Cluster_Labels[species == "Stephensi", 1]
Cluster_Identity_Sub <- Cluster_Labels[species == "Subpictus", 1]

STAN_Data <- list(Number_Clusters = Number_Clusters, Number_Envt_Variables = Number_Envt_Variables, 
                  Number_Obs_Ann = Number_Obs_Ann, Number_Obs_Cul = Number_Obs_Cul, Number_Obs_Dir = Number_Obs_Dir,
                  Number_Obs_Fluv = Number_Obs_Fluv, Number_Obs_Min = Number_Obs_Min, Number_Obs_Ste = Number_Obs_Ste, Number_Obs_Sub = Number_Obs_Sub,
                  Envt_Variables_Ann = Envt_Variables_Ann, Envt_Variables_Cul = Envt_Variables_Cul, Envt_Variables_Dir = Envt_Variables_Dir,
                  Envt_Variables_Fluv = Envt_Variables_Fluv, Envt_Variables_Min = Envt_Variables_Min, Envt_Variables_Ste = Envt_Variables_Ste, Envt_Variables_Sub = Envt_Variables_Sub,
                  Cluster_Identity_Ann = Cluster_Identity_Ann, Cluster_Identity_Cul = Cluster_Identity_Cul, Cluster_Identity_Dir = Cluster_Identity_Dir,
                  Cluster_Identity_Fluv = Cluster_Identity_Fluv, Cluster_Identity_Min = Cluster_Identity_Min, Cluster_Identity_Ste = Cluster_Identity_Ste, Cluster_Identity_Sub = Cluster_Identity_Sub)

STAN_GLM <- stan_model("Model_Files/Multinomial_Spec_Spec_Intercepts.stan")
options(mc.cores = 3)
STAN_fit <- sampling(STAN_GLM, data = STAN_Data, iter = 3000, chains = 3, refresh = 500, control = list(max_treedepth = 15, adapt_delta = 0.95))

check_hmc_diagnostics(STAN_fit)
fit <- rstan::extract(STAN_fit)

if (prior == "informative") {
  saveRDS(fit, "Outputs/Logistic_Regression_Output/Informative_Prior/Reduced_Subset_STAN_Output.rds")
} else if (prior == "uninformative") {
  saveRDS(fit, "Outputs/Logistic_Regression_Output/Uninformative_Prior/Reduced_Subset_STAN_Output.rds")
}
fit <- readRDS("Outputs/Logistic_Regression_Output/Informative_Prior/Reduced_Subset_STAN_Output.rds")


#######################################################################################################
##                                                                                                   ##
##                         Analysing and Exploring Logistic Regression Results                       ## 
##                                                                                                   ##
#######################################################################################################
spec_coefs <- rbind(apply(fit$alpha_ann, 2, mean), apply(fit$alpha_cul, 2, mean), apply(fit$alpha_dir, 2, mean), 
                    apply(fit$alpha_fluv, 2, mean), apply(fit$alpha_min, 2, mean), apply(fit$alpha_ste, 2, mean), apply(fit$alpha_sub, 2, mean))
row.names(spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
cluster_spec_melt <- melt(spec_coefs)
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = rev(unique(cluster_spec_melt$Var1[order(cluster_spec_melt$Var1)])))
ggplot(cluster_spec_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")

envt_betas <- matrix(nrow = number_covariates, ncol = Number_Clusters)
for (i in 1:number_covariates) {
  envt_betas[i, ] <- apply(fit$envt_betas[ , i, ], 2, mean)
}
row.names(envt_betas) <- variable_names

envt_melt_test <- melt(envt_betas)
base_size <- 9
ggplot(envt_melt_test, aes(x = Var2, y = Var1)) + 
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
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme_grey()


#######################################################################################################
##                                                                                                   ##
##              Assessing Predictive Performance of Logistic Regression Output                       ## 
##                                                                                                   ##
#######################################################################################################
pred_prob_Ann <- predict_probabilities(Cluster_Identity_Ann, "Cluster_Probabilities_Ann", fit, Number_Clusters)
pred_prob_Cul <- predict_probabilities(Cluster_Identity_Cul, "Cluster_Probabilities_Cul", fit, Number_Clusters)
pred_prob_Dir <- predict_probabilities(Cluster_Identity_Dir, "Cluster_Probabilities_Dir", fit, Number_Clusters)
pred_prob_Fluv <- predict_probabilities(Cluster_Identity_Fluv, "Cluster_Probabilities_Fluv", fit, Number_Clusters)
pred_prob_Min <- predict_probabilities(Cluster_Identity_Min, "Cluster_Probabilities_Min", fit, Number_Clusters)
pred_prob_Ste <- predict_probabilities(Cluster_Identity_Ste, "Cluster_Probabilities_Ste", fit, Number_Clusters)
pred_prob_Sub <- predict_probabilities(Cluster_Identity_Sub, "Cluster_Probabilities_Sub", fit, Number_Clusters)

pred_clus_Ann <- predict_cluster(pred_prob_Ann)
pred_clus_Cul <- predict_cluster(pred_prob_Cul)
pred_clus_Dir <- predict_cluster(pred_prob_Dir)
pred_clus_Fluv <- predict_cluster(pred_prob_Fluv)
pred_clus_Min <- predict_cluster(pred_prob_Min)
pred_clus_Ste <- predict_cluster(pred_prob_Ste)
pred_clus_Sub <- predict_cluster(pred_prob_Sub)

Correct_Ann <- calculate_correct_predictions(Cluster_Identity_Ann, pred_clus_Ann)
Correct_Cul <- calculate_correct_predictions(Cluster_Identity_Cul, pred_clus_Cul)
Correct_Dir <- calculate_correct_predictions(Cluster_Identity_Dir, pred_clus_Dir)
Correct_Fluv <- calculate_correct_predictions(Cluster_Identity_Fluv, pred_clus_Fluv)
Correct_Min <- calculate_correct_predictions(Cluster_Identity_Min, pred_clus_Min)
Correct_Ste <- calculate_correct_predictions(Cluster_Identity_Ste, pred_clus_Ste)
Correct_Sub <- calculate_correct_predictions(Cluster_Identity_Sub, pred_clus_Sub)

total_correct <- Correct_Ann + Correct_Cul + Correct_Dir + Correct_Fluv + Correct_Min + Correct_Ste + Correct_Sub
correct_proportion <- total_correct/272

#######################################################################################################
##                                                                                                   ##
##                    Analysing and Exploring Logistic Regression Results for SPECIES                ## 
##                                                                                                   ##
#######################################################################################################
spec_coefs <- rbind(apply(fit$alpha_ann, 2, mean), apply(fit$alpha_cul, 2, mean), apply(fit$alpha_dir, 2, mean), 
                    apply(fit$alpha_fluv, 2, mean), apply(fit$alpha_min, 2, mean), apply(fit$alpha_ste, 2, mean), apply(fit$alpha_sub, 2, mean))
row.names(spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
colnames(spec_coefs) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

species_coefficient_clustering <- hclust(dist(spec_coefs))
spec_dend <- species_coefficient_clustering %>%
  as.dendrogram() %>%
  color_branches(k = 7, col = c("#898989", "#2EBEEA",  "#A54D2C", "#29B200" , "#E0521A", "#A948EA", "#E8A50B")) %>%
  color_labels(col = c("#898989", "#2EBEEA",  "#A54D2C", "#29B200" , "#E0521A", "#A948EA", "#E8A50B")) %>%
  set("labels_cex", 1.5) %>%
  set("highlight_branches_lwd", 3)
plot(spec_dend, las = 1, ylab = "", xlab = "", main = "", horiz = TRUE, lwd = 3, cex = 2)

cluster_spec_melt <- melt(spec_coefs)
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = rev(unique(cluster_spec_melt$Var1[order(cluster_spec_melt$Var1)])))
ggplot(cluster_spec_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red", 
                       limits = c(-1, 1.5), breaks = c(-0.75, 0, 0.75, 1.5)) + 
  xlab ("") + 
  ylab("") +
  theme_minimal() +
  theme(text = element_text(size = 20))

spec_melt_cor <- melt(cor(spec_coefs))
spec_melt_cor$Var1 <- factor(spec_melt_cor$Var1, levels = rev(unique(spec_melt_cor$Var1[order(spec_melt_cor$Var1)])))
base_size <- 9
ggplot(spec_melt_cor, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")


#######################################################################################################
##                                                                                                   ##
##              Analysing and Exploring Logistic Regression Results for ENVIRONMENTAL                ## 
##                                                                                                   ##
#######################################################################################################
envt_betas <- matrix(nrow = length(variable_names), ncol = 4)
for (i in 1:length(variable_names)) {
  envt_betas[i, ] <- apply(fit$envt_betas[ , i, ], 2, mean)
}
row.names(envt_betas) <- variable_names
colnames(envt_betas) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

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
upset(fromList(listInput), keep.order = TRUE, intersections = intersections, set_size.show = FALSE,
      order.by = "freq", point.size = 7.5, line.size = 3.5, 
      mainbar.y.label = "Top Ecological Variable Intersections",
      text.scale = 1.6)

# Supplementary Figure - Ecological Coefficient Value Heatmap
envt_melt <- melt(envt_betas)
envt_melt$Var1 <- factor(envt_melt$Var1, levels = rev(unique(envt_melt$Var1[order(envt_melt$Var1)])))
base_size <- 9
ggplot(envt_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red") +
  xlab ("") + 
  ylab("") +
  theme_minimal()


# Supplementary Figure - Ecological Coefficient Correlation Across Clusters Plot 
envt_melt_cor <- melt(cor(envt_betas))
envt_melt_cor$Var1 <- factor(envt_melt_cor$Var1, levels = rev(unique(envt_melt_cor$Var1[order(envt_melt_cor$Var1)])))
base_size <- 9
ggplot(envt_melt_cor, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")



#####################
### Miscellaneous ###
#####################
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = unique(cluster_spec_melt$Var1))
ggplot(cluster_spec_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  xlab ("") + 
  ylab("Coefficient Size") +
  theme_minimal() +
  theme(text = element_text(size=20), axis.title.y = element_text(vjust = 5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  guides(fill=guide_legend(title="Species")) +
  scale_fill_manual(values = c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C"))

# Supplementary Figure - Cluster 1 Envt Coefs Ordered
cluster_1_betas <- envt_betas[, 1]
index_for_ordering <- order(cluster_1_betas)
envt_melt_test <- melt(envt_betas)
envt_melt_test$Var1 <- factor(envt_melt_test$Var1, levels = envt_melt_test$Var1[index_for_ordering])
base_size <- 9
ggplot(envt_melt_test, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("")


envt_betas[, 1][order(abs(envt_betas[, 1]), decreasing = TRUE)[1:5]]
envt_betas[, 2][order(abs(envt_betas[, 2]), decreasing = TRUE)[1:5]]
envt_betas[, 3][order(abs(envt_betas[, 3]), decreasing = TRUE)[1:5]]
envt_betas[, 4][order(abs(envt_betas[, 4]), decreasing = TRUE)[1:5]]
top_vars <- c("Water Areas Max Extent", "LC8", "EVI Mean", "Mean Temp Wettest Quarter", "Rain Warmest Quarter", 
              "DCW Distance to Water", "PET Yearly Average", "Flow Accumulation", "Mean Temp Coldest Quarter", 
              "Night LST SD", "WC A3", "Urban Footprint", "Rain Coldest Quarter", "Tasseled Cap Brightness Mean", 
              "Rain Coldest Quarter", "India Pop Density 2010", "WC P3", "Max Temp Warmest Month")
envt_melt_top <- envt_melt[envt_melt$Var1 %in% top_vars, ] 
envt_melt_top$Var1 <- factor(envt_melt_top$Var1, levels = rev(unique(envt_melt_top$Var1[order(envt_melt_top$Var1)])))
base_size <- 9
ggplot(envt_melt_top, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red") +
  xlab ("") + 
  ylab("") +
  theme_minimal()



