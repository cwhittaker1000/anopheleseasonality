#######################################################################################################
##                                                                                                   ##
##                         Initial Loading of Libraries, Shapefiles & Data                           ##
##                                                                                                   ##
#######################################################################################################
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/")
source("R_Scripts/1_Covariate_Extraction_and_Collation/2_Polygon_Construction.R")
source("Functions/Periodic_Kernel_GP_Fitting_Functions.R")
source("Functions/Von_Mises_Fitting_Functions.R")
source("Functions/Time_Series_Operation_Functions.R")
adm1data <- getData('GADM', country = 'IND', level = 1)
mosquito_data <- read.csv("Datasets/Systematic_Review/Processed_Catch_Data.csv", stringsAsFactors = FALSE)
keep <- mosquito_data$Keep
species <- mosquito_data$Species
time_series <- mosquito_data

#######################################################################################################
##                                                                                                   ##
##                        Plotting the Collated Time Series On a Map of India                        ## 
##                                                                                                   ##
#######################################################################################################
par(mfrow = c(1, 1))
species_geos <- data.frame(x = NA, y = NA, species = NA, stringsAsFactors = FALSE)
counter <- 1
for (i in 1:272) {
  
  index <- time_series$Location_ID[i]
  matching_location_index <- which(location_coordinates$Location_ID %in% index)
  
  if (length(matching_location_index) == 1) {
    coords <- st_coordinates(location_coordinates$geometry[matching_location_index])[, 1:2]
    if (length(coords) == 2) {
      centroid <- coords
    } else {
      centroid <- geosphere::centroid(coords)
    }
    species_geos[counter, 1] <- centroid[1]
    species_geos[counter, 2] <- centroid[2]
    species_geos[counter, 3] <- time_series$Species[i]
    counter <- counter + 1
    
  } else {
    print(length(matching_location_index))
    for (j in 1:length(matching_location_index)) {
      coords <- st_coordinates(location_coordinates$geometry[matching_location_index[j]])[, 1:2]
      if (length(coords) == 2) {
        centroid <- coords
      } else {
        centroid <- geosphere::centroid(coords)
      }
      species_geos[counter, 1] <- centroid[1]
      species_geos[counter, 2] <- centroid[2]
      species_geos[counter, 3] <- time_series$Species[i]
      counter <- counter + 1
    }
  }

}

set.seed(44)
palette(c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C"))
plot(adm1data, col = "#FFF6E2", lwd = 2)
points(species_geos[, 1] + rnorm(length(species_geos[, 1]), 0, 0.2), species_geos[, 2] + rnorm(length(species_geos[, 2]), 0, 0.2), col = as.factor(species_geos[, 3]), pch = 20, cex = 2)
pdf("Figures/Figure_1/Figure_1A_Indian_Map.pdf", height = 8, width = 10)
palette(c("#29B200" , "#2EBEEA", "#A948EA", "#898989", "#E8A50B", "#E0521A", "#A54D2C"))
plot(adm1data, col = "#FFF6E2", lwd = 2)
points(species_geos[, 1] + rnorm(length(species_geos[, 1]), 0, 0.2), species_geos[, 2] + rnorm(length(species_geos[, 2]), 0, 0.2), col = as.factor(species_geos[, 3]), pch = 20, cex = 2)
dev.off()

