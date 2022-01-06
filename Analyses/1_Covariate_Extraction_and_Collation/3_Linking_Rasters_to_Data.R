#######################################################################################################
##                                                                                                   ##
##                Loading Required Libraries, Raster Data & Setting the Working Directory            ##
##                                                                                                   ##
#######################################################################################################
library(rasterVis); library(raster); library(sp); library(classInt); library(colorspace);
library(maps); library(sp); library(rgeos); library(rgdal); library(maptools); library(dplyr);
library(tidyr); library(tmap); library(maps); library(scales); library(mapproj); library(gdalUtils);
library(latticeExtra); library(geosphere); library(sf); library(ggplot2); library(here)
library(tsne); library(FactoMineR); library(factoextra); library(MALDIquant); library(rgl)

source("Analyses/1_Covariate_Extraction_and_Collation/1_Loading_Rasters.R")
source("Analyses/1_Covariate_Extraction_and_Collation/2_Polygon_Construction.R")
India_Surroundings_Extent <- extent(60, 100, 0, 40) # Creates extent covering India - use to crop raster


########################################################################################################
##                                                                                                    ##
##                       Stacking All the Rasters Together For Value Extraction                       ##
##                                                                                                    ##
##           Projecting all rasters to the same resolution and stacking them together.                ##
##                                                                                                    ##
########################################################################################################
rasters <- c("Annual_Mean_Temperature", "Mean_Diurnal_Range", "Isothermality", "Temperature_Seasonality", 
             "Max_Temp_Warmest_Month", "Min_Temp_Coldest_Month", "Temp_Annual_Range", "Mean_Temp_Wettest_Quarter", 
             "Mean_Temp_Driest_Quarter", "Mean_Temp_Warmest_Quarter", "Mean_Temp_Coldest_Quarter", "Annual_Rain", 
             "Rain_Wettest_Month", "Rain_Driest_Month", "Rain_Seasonality", "Rain_Wettest_Quarter", "Rain_Driest_Quarter", 
             "Rain_Warmest_Quarter", "Rain_Coldest_Quarter",            
             "PET_Yearly_Average", "Aridity_Yearly_Average",
             "India_Pop_Density_2010",
             "Day_LST_Mean", "Day_LST_SD", "Night_LST_Mean", "Night_LST_SD",                                           
             "Tasseled_Cap_Wetness_Mean", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_Mean", "Tasseled_Cap_Brightness_SD",
             "Elevation", "Specific_Humidity_Mean", "Specific_Humidity_SD", "EVI_Mean", "Flow_Accumulation", 
             "Water_Areas_Max_Extent", "Water_Areas_Seasonality", "Water_Areas_Occurrence", "Water_Areas_Recurrence",
             "DCW_Distance_to_Water", "WWF_Distance_to_Water", 
             "City_Accessibility", 
             "CHIRPS_Max", "CHIRPS_Min", "CHIRPS_Mean",                                                                                                                                                     
             "WC_A0", "WC_A1", "WC_A2", "WC_A3", "WC_P0", "WC_P1", "WC_P2", "WC_P3", 
             "Urban_Footprint", "Irrigated_Areas", "Dominant_Landcover") #"GlobCover", "Globcover_Wet_Areas", 

number_of_covariates <- length(rasters)
raster_stack <- stack()
for (i in 1:number_of_covariates) { 
  if (i != number_of_covariates) {
    projected_raster <- projectRaster(eval(parse(text = rasters[i])), Day_LST_Mean)
    cropped_raster_layer <- crop(projected_raster, India_Surroundings_Extent) 
    cropped_raster_layer <- scale(cropped_raster_layer) 
    raster_stack <- stack(raster_stack, cropped_raster_layer)
  } else {
    projected_raster <- projectRaster(eval(parse(text = rasters[i])), Day_LST_Mean, method = "ngb") # required for categorical landcover
    cropped_raster_layer <- crop(projected_raster, India_Surroundings_Extent) 
    raster_stack <- stack(raster_stack, cropped_raster_layer)
  }
  print(i)
}
names(raster_stack) <- rasters
writeRaster(raster_stack,"Datasets/Extracted_Covariates_for_Modelling/Raster_Stack.grd", format = "raster", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##       Extracting Raster Values for Each Locations and Normalising (Mean = 0, Unit Variance)       ##
##                                                                                                   ##
##    Extracting raster values for each of the locations where mosquitoes were sampled. The          ##
##    following is carried out:                                                                      ##
##                                                                                                   ## 
##        1) Manual processing/dealing with instances where polygons are slightly over               ##
##           water and hence would return NA values in many instances when we average.               ##
##        2) Extraction and collation of raster values for all of the different locations.           ##
##        3) Manual processing for locations comprised of 2 distinct areas In these cases,           ##
##           the extracted covariate value is averaged over the 2 locations.                         ##
##                                                                                                   ##
#######################################################################################################
# Manual Processing of Location 23 to Move it Marginally More Inland
location_coordinates$geometry[23] <- st_geometry(st_point(c(80.2600, 13.0002))) 

# Extracting Covariate Values for Both Points and Locations With >1 Associated Pair of Geocoordinates
raw_covariate_values <- matrix(nrow = length(location_coordinates$Location_ID), ncol = number_of_covariates)
for (i in 1:length(location_coordinates$geometry)) {
  for (j in 1:number_of_covariates) {
    if (j == number_of_covariates) {  # need to process landcover, which is the last covariate, differently
      if (length(st_coordinates(location_coordinates$geometry[i])[, 1]) == 1) {
        raw_covariate_values[i, j] <- raster::extract(raster_stack[[j]], st_coordinates(location_coordinates$geometry[i]))
      }
      else {  # Calculate all the raster values that a polygon spans then calculates the mean
        coordinates_for_the_polygon <- st_coordinates(location_coordinates$geometry[i])[, 1:2]
        temp_polygon <- spPolygons(coordinates_for_the_polygon)
        polygon_covariate_values <- unlist(raster::extract(raster_stack[[j]], temp_polygon))
        raw_covariate_values[i, j] <- get_mode(na.omit(polygon_covariate_values))
      } 
    }
    else {
      if (length(st_coordinates(location_coordinates$geometry[i])[, 1]) == 1) {
        raw_covariate_values[i, j] <- raster::extract(raster_stack[[j]], st_coordinates(location_coordinates$geometry[i]))
      }
      else {  # Calculate all the raster values that a polygon spans then calculates the mean
        coordinates_for_the_polygon <- st_coordinates(location_coordinates$geometry[i])[, 1:2]
        temp_polygon <- spPolygons(coordinates_for_the_polygon)
        polygon_covariate_values <- unlist(raster::extract(raster_stack[[j]], temp_polygon))
        raw_covariate_values[i, j] <- mean(na.omit(polygon_covariate_values), na.rm = TRUE)
      } 
    }
  }
  print(i)
}
colnames(raw_covariate_values) <- rasters

# Checking for NAs and Which Columns They're In
for (i in 1:dim(raw_covariate_values)[2]) {
  NAs <- sum(is.na(raw_covariate_values[, i]))
  print(paste0("Column ", i, " has ", NAs, " NAs" ))
}

# Replacing the two sets of NAs (for only two locations) in Specific Humidity Mean + SD with 0
raw_covariate_values[25:26, 32:33] <- 0

# Taking the Mean of the Covariate Values for For Each Polyon That Has 2 Points Comprising It
# Then Removing Every Other Row (These Rows Are Now Redundant) - this leaves 1 row per pair of points,
# with covariate values that are the mean of the covariate values for each of the pairs
two_point_indices <- which(location_coordinates$Geo_Type == "2_Singles") 
counter_index <- min(two_point_indices) # starts the counter at where the polygons with two geoocordinate pairs begins in the dataframe
for (i in two_point_indices) {
  if (i %% 2 == 0) { # selects every other row (specifically, the first row of each pair)
    for (j in 1:(number_of_covariates)) {
      raw_covariate_values[counter_index, j] <- mean(c(raw_covariate_values[i, j], raw_covariate_values[i + 1, j]))
    }
    counter_index <- counter_index + 1
  }
}
covariate_rows_to_be_removed <- 107:117
covariate_values <- raw_covariate_values[-covariate_rows_to_be_removed, ]
locations_IDs_to_be_removed <- seq(min(two_point_indices) + 1, max(two_point_indices), 2)
location_IDs <- location_coordinates$Location_ID[-locations_IDs_to_be_removed] # note different rows removed - this is because
                                                                               # when the raw_covariate_values pairs of locations
                                                                               # are processed, we sequentially replace the next
                                                                               # row in the dataset, see Lines 102-111 to see 

covariate_values <- cbind(location_IDs, covariate_values)
saveRDS(covariate_values, file = "Datasets/Extracted_Covariates_for_Modelling/Extracted_Covariates.rds")

