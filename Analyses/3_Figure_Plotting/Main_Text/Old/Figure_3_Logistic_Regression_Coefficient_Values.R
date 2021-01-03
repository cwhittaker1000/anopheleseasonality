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
  fit <- readRDS("Outputs/Logistic_Regression_Output/Informative_Prior/STAN_Output.rds")
} else if (prior == "uninformative") {
  fit <- readRDS("Outputs/Logistic_Regression_Output/Uninformative_Prior/STAN_Output.rds")
}

#######################################################################################################
##                                                                                                   ##
##                    Analysing and Exploring Logistic Regression Results for SPECIES                ## 
##                                                                                                   ##
#######################################################################################################
spec_coefs <- rbind(apply(fit$alpha_ann, 2, mean), apply(fit$alpha_cul, 2, mean), apply(fit$alpha_dir, 2, mean), 
                    apply(fit$alpha_fluv, 2, mean), apply(fit$alpha_min, 2, mean), apply(fit$alpha_ste, 2, mean), apply(fit$alpha_sub, 2, mean))
row.names(spec_coefs) <- c("Annularis", "Culicifacies", "Dirus", "Fluviatilis", "Minimus", "Stephensi", "Subpictus")
colnames(spec_coefs) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

# Figure 3A - Hierachical Clustering of Species Coefficients
species_coefficient_clustering <- hclust(dist(spec_coefs))
spec_dend <- species_coefficient_clustering %>%
  as.dendrogram() %>%
  color_branches(k = 7, col = c("#898989", "#2EBEEA",  "#A54D2C", "#29B200" , "#E0521A", "#A948EA", "#E8A50B")) %>%
  color_labels(col = c("#898989", "#2EBEEA",  "#A54D2C", "#29B200" , "#E0521A", "#A948EA", "#E8A50B")) %>%
  set("labels_cex", 1.5) %>%
  set("highlight_branches_lwd", 3)
pdf(file = "Figures/Figure_3A_Species_Coefficient_Dendogram.pdf", width = 6.5, height = 5.5)
plot(spec_dend, las = 1, ylab = "", xlab = "", main = "", horiz = TRUE, lwd = 3, cex = 2)
dev.off()

# Figure 3B - Species Coefficient Values for Each Cluster
cluster_spec_melt <- melt(spec_coefs)
cluster_spec_melt$Var1 <- factor(cluster_spec_melt$Var1, levels = rev(unique(cluster_spec_melt$Var1[order(cluster_spec_melt$Var1)])))
pdf(file = "Figures/Figure_3B_Species_Coefficient_Heatmap.pdf", width = 6.5, height = 5.5)
ggplot(cluster_spec_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2("Effect Size", low = "blue", mid = "white", high = "red", 
                       limits = c(-0.55, 1), breaks = c(-0.5, 0, 0.5, 1)) + 
  xlab ("") + 
  ylab("") +
  theme_minimal() +
  theme(text = element_text(size = 20))
dev.off()

# Supplementary Figure - Species Coefficient Correlations Across Clusters 
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
variable_names <- c("Annual Mean Temperature", "Mean Diurnal Range", "Isothermality", "Temperature Seasonality", 
                    "Max Temp Warmest Month", "Min Temp Coldest Month", "Temp Annual Range", "Mean Temp Wettest Quarter", 
                    "Mean Temp Driest Quartest", "Mean Temp Warmest Quarter", "Mean Temp Coldest Quarter", "Annual Rain", 
                    "Rain Wettest Month", "Rain Driest Month", "Rain Seasonality", "Rain Wettest Quarter", "Rain Driest Quarter", 
                    "Rain Warmest Quarter", "Rain Coldest Quarter",            
                    "PET Yearly Average", "Aridity Yearly Average",
                    "India Pop Density 2010",
                    "Day LST Mean", "Day LST SD", "Night LST Mean", "Night LST SD",                                           
                    "Tasseled Cap Wetness Mean", "Tasseled Cap Wetness SD", "Tasseled Cap Brightness Mean", "Tasseled Cap Brightness SD",
                    "Elevation", "Specific Humidity Mean", "Specific Humidity SD", "EVI Mean", "Flow Accumulation", 
                    "Water Areas Max Extent", "Water Areas Seasonality", "Water Areas Occurrence", "Water Areas Recurrence",
                    "DCW Distance to Water", "WWF Distance to Water", 
                    "City Accessibility", 
                    "CHIRPS Max", "CHIRPS Min", "CHIRPS Mean",                                                                                                                                                     
                    "WC A0", "WC A1", "WC A2", "WC A3", "WC P0", "WC P1", "WC P2", "WC P3", 
                    "Urban Footprint", "Irrigated Areas", "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8", "LC9", "LC10", "LC11")

envt_betas <- matrix(nrow = 66, ncol = 4)
for (i in 1:66) {
  envt_betas[i, ] <- apply(fit$envt_betas[ , i, ], 2, mean)
}
row.names(envt_betas) <- variable_names
colnames(envt_betas) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

# Figure 3C - UpsetR Plot of Intersections of Top Ecological Coefficients for Each Cluster Pair 
one_top_five <- envt_betas[, 1][order(abs(envt_betas[, 1]), decreasing = TRUE)][1:5]
two_top_five <- envt_betas[, 2][order(abs(envt_betas[, 2]), decreasing = TRUE)][1:5]
three_top_five <- envt_betas[, 3][order(abs(envt_betas[, 3]), decreasing = TRUE)][1:5]
four_top_five <- envt_betas[, 4][order(abs(envt_betas[, 4]), decreasing = TRUE)][1:5]


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

