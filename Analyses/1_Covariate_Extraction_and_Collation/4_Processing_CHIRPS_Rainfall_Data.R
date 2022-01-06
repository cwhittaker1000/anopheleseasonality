#######################################################################################################
##                                                                                                   ##
##      Loading Required Libraries, Functions, Scripts & Setting the Working Directory               ##
##                                                                                                   ##
##    The code below sources code that processes the location data for each mosquito sampling        ##
##    site, and loads in other information about the timing and duration of the sampling period.     ##
##    Together, this information is used to extract daily CHIRPS rainfall data for a given           ##
##    location over the entirety of the sampling period. This rainfall data is stored in .csv        ##
##    files, with a single .csv file for each sampling site.                                         ##
##                                                                                                   ##
#######################################################################################################
library(dplyr); library(raster); library(rgdal); library(sf); library(raster); library(sp); library(tictoc); library(here)

source("Functions/CHIRPS_Rainfall_Processing_Functions.R") # functions
source("Analyses/1_Covariate_Extraction_and_Collation/2_Polygon_Construction.R") # spatial information
spatial_information <- read.csv("Datasets/Systematic_Review/Geographical_Data.csv", stringsAsFactors = FALSE) # temporal information

spatial_metadata <- spatial_information[, c(3, 7, 8)]
location_coordinates <- location_coordinates[order(location_coordinates$Location_ID), ]
directory <- paste0(here("Datasets/CHIRPS_Rainfall_Data/"), "/Overall_India_CHIRPS_Rainfall")

# Extract Location and Date Specific Rainfall for All Locations 
for (i in 1:max(location_coordinates$Location_ID)) {
  
  # Creating the Required Variables to Extract the Data (Extract With Year Buffer Either Side of Sampling Period)
  years <- seq(spatial_metadata$Sampling.Start[i] - 1, spatial_metadata$Sampling.End[i] + 1, 1)
  location_ID <- i

  # Extracts the Rainfall Time Series
  raw_rainfall_time_series_output <- rain_time_series_generator(years, location_ID, location_coordinates, directory) 
  rainfall_time_series <- raw_rainfall_time_series_output$rainfall_time_series
  dates_for_csv <- raw_rainfall_time_series_output$dates
  
  # Creates Time Series Output and Saves as .csv
  time_series_output <- as.data.frame(dates_for_csv)
  time_series_output[, 2] <- rainfall_time_series
  colnames(time_series_output) <- c("Date", "Rainfall") 
  
  write.csv(time_series_output, file = paste0("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Indian_Vectors/1. Temporal Vector Dynamics/Data/Location_Specific_CHIRPS_Rainfall/LocID_", 
                                              location_metadata$Location_ID[i], "_RefID_", location_metadata$Ref_ID[i], ".csv", sep = ""), row.names = FALSE)
  
  print(paste0("Location ID ", i, " has been completed", sep = ""))
  
}

