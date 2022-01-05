#######################################################################################################
##                                                                                                   ##
##                          Loading Required Libraries and Other Variables                           ##
##                                                                                                   ##
#######################################################################################################
library(rasterVis); library(raster); library(sp); library(classInt); library(colorspace); library(maps); 
library(sp); library(rgeos); library(rgdal); library(maptools); library(dplyr); library(tidyr); 
library(tmap); library(maps); library(scales); library(mapproj); library(gdalUtils); library(latticeExtra); 
library(geosphere); library(sf); library(viridis); library(gdistance); library(tictoc)
source("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/Functions/Spatial_Data_Processing_Functions.R")
source("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/Functions/CHIRPS_Rainfall_Processing_Functions.R")
India_Surroundings_Extent <- extent(60, 100, 0, 40) # Creates extent covering India - use to crop raster
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/Datasets/Satellite_Covariates/")


#######################################################################################################
##                                                                                                   ##
##                                       Bioclimatic Variables                                       ##
##                                                                                                   ##
##   Bioclimatic variables derived from monthly temperature and rainfall values in order to          ##
##   generate variables that are biologically relevant/meaningful (e.g. Annual Mean Temperature,     ##
##   Mean Diurnal Range etc). Based on the WorldClim data.                                           ##
##     Note: Some variables might not require standardisation as the others do.                      ##
##     Available Here: http://worldclim.org/version2                                                 ##
##                                                                                                   ##
#######################################################################################################
lat_vec <- c(30, 60, 30, 60) # tiles for Bioclimatic Variables - Latitude Coordinates
lon_vec <- c(60, 60, 90, 90) # tiles for Bioclimatic Variables - Longitude Coordinates
bioclim_raster_1 <- getData('worldclim', var = 'bio', res = 0.5, lon = lon_vec[1], lat = lat_vec[1], path = "Raw/")
bioclim_raster_2 <- getData('worldclim', var = 'bio', res = 0.5, lon = lon_vec[2], lat = lat_vec[2], path = "Raw/")
bioclim_raster_3 <- getData('worldclim', var = 'bio', res = 0.5, lon = lon_vec[3], lat = lat_vec[3], path = "Raw/")
bioclim_raster_4 <- getData('worldclim', var ='bio', res = 0.5, lon = lon_vec[4], lat = lat_vec[4], path = "Raw/")
BioClimNames <- c("Annual_Mean_Temperature", "Mean_Diurnal_Range", "Isothermality", 
                  "Temperature_Seasonality", "Max_Temp_Warmest_Month", "Min_Temp_Coldest_Month", 
                  "Temp_Annual_Range", "Mean_Temp_Wettest_Quarter", "Mean_Temp_Driest_Quarter", 
                  "Mean_Temp_Warmest_Quarter", "Mean_Temp_Coldest_Quarter", "Annual_Rain", 
                  "Rain_Wettest_Month", "Rain_Driest_Month", "Rain_Seasonality", 
                  "Rain_Wettest_Quarter", "Rain_Driest_Quarter", "Rain_Warmest_Quarter", "Rain_Coldest_Quarter")
for (i in 1:19) {
  x <- paste(BioClimNames[i])
  mosaic_variable <- mosaic(bioclim_raster_1[[i]], bioclim_raster_2[[i]], bioclim_raster_3[[i]], bioclim_raster_4[[i]], fun = mean)
  eval(call("<-", as.name(x), mosaic_variable))
  print(i)
}
raster::writeRaster(Annual_Mean_Temperature, filename = "Processed/Annual_Mean_Temperature.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Mean_Diurnal_Range, filename = "Processed/Mean_Diurnal_Range.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Isothermality, filename = "Processed/Isothermality.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Temperature_Seasonality, filename = "Processed/Temperature_Seasonality.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Max_Temp_Warmest_Month, filename = "Processed/Max_Temp_Warmest_Month.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Min_Temp_Coldest_Month, filename = "Processed/Min_Temp_Coldest_Month.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Temp_Annual_Range, filename = "Processed/Temp_Annual_Range.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Mean_Temp_Wettest_Quarter, filename = "Processed/Mean_Temp_Wettest_Quarter.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Mean_Temp_Driest_Quartest, filename = "Processed/Mean_Temp_Driest_Quartest.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Mean_Temp_Warmest_Quarter, filename = "Processed/Mean_Temp_Warmest_Quarter.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Mean_Temp_Coldest_Quarter, filename = "Processed/Mean_Temp_Coldest_Quarter.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Annual_Rain, filename = "Processed/Annual_Rain.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Wettest_Month, filename = "Processed/Rain_Wettest_Month.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Driest_Month, filename = "Processed/Rain_Driest_Month.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Seasonality, filename = "Processed/Rain_Seasonality.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Wettest_Quarter, filename = "Processed/Rain_Wettest_Quarter.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Driest_Quarter, filename = "Processed/Rain_Driest_Quarter.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Warmest_Quarter, filename = "Processed/Rain_Warmest_Quarter.tif", format = "GTiff", overwrite = TRUE)
raster::writeRaster(Rain_Coldest_Quarter, filename = "Processed/Rain_Coldest_Quarter.tif", format = "GTiff", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##                               Potential Evapotranspiration & Aridity                              ##
##                                                                                                   ##
##  Global Potential Evapotranspiration (Annual Average 1950-2000)                                   ##
##    Potential Evapotranspiration (PET) is a measure of the ability of the atmosphere to remove     ##
##    water through the ET process. I.e. describes the amount of evaporation that would occur if     ##
##    a sufficient water source were available. PET measures dryness to an extent. Based on          ## 
##    WorldClim data.                                                                                ##
##       Available Here: https://cgiarcsi.community/data/global-aridity-and-pet-database/            ##
##                                                                                                   ##
##  Global Aridity Index (Annual Average 1950-2000)                                                  ##
##    Measures the aridity of an area, which broadly corresponds to how extenstively it lack         ##
##    the necessary moisture to achieve positive vegetation growth. How arid an area is.             ##
##    Based on WorldClim data.                                                                       ##
##      Note: Readme says need to divide by 10,000 to get the values in the correct units.           ##
##      Available Here: https://cgiarcsi.community/data/global-aridity-and-pet-database/             ##
##                                                                                                   ##
#######################################################################################################
PET_Yearly_Average_Raw <- raster("Raw/PET_Annual/pet_he_yr/w001001.adf")
PET_Yearly_Average <- crop(PET_Yearly_Average_Raw, India_Surroundings_Extent)
raster::writeRaster(PET_Yearly_Average, filename = "Processed/PET_Annual.tif", format = "GTiff", overwrite = TRUE)
Aridity_Yearly_Average_Raw <- raster("Raw/AI_Annual/ai_yr/w001001.adf")
Aridity_Yearly_Average <- crop(Aridity_Yearly_Average_Raw, India_Surroundings_Extent)
Aridity_Yearly_Average <- Aridity_Yearly_Average * 0.0001
raster::writeRaster(Aridity_Yearly_Average, filename = "Processed/GAI_Annual.tif", format = "GTiff", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##                                      Population Density                                           ##
##   Population Density (2010 currently)                                                             ##
##      Raster maps detailing the population per pixel from WorldPop's website. WorldPop             ##
##      offers options for either 2010, 2011, 2015 or 2020 (projected), and either as their raw      ##
##      estimates or as estimates once national total have been adjusted to match UN population      ##
##      division estimates.                                                                          ##
##          Note: Unclear how useful this will be given the large timespan my data lies over.        ##
##          Available from: http://www.worldpop.org.uk/data/files/                                   ##
##                                                                                                   ##
#######################################################################################################
India_Pop_Density_2010 <-  raster("Raw/IND_ppp_2010_adj_v2.tif")
India_Pop_Density_2010 <- aggregate(x = India_Pop_Density_2010, fact = 10, fun = sum)
India_Pop_Density_2010 <- projectRaster(India_Pop_Density_2010, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ')
raster::writeRaster(India_Pop_Density_2010, filename = "Processed/Pop_Density_Aggregated_Tenfold.tif", format = "GTiff", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##                     Enhanced Vegetation Index (Yearly Average 2001-2015)                          ## 
##                                                                                                   ##
##                                                                                                   ##
##      Enhanced Vegetation Index (Yearly Average 2001-2015)                                         ## 
##          Raster maps detailed the enhanced vegetation index value per pixel from Google Earth     ##
##          Engine's website. Originally derived by the Malaria Atlas Project and added to GEE as    ##
##          a dataset that can be manipulated. Further details on the website. Yearly averages       ##
##          for the years 2001 through to 2015.                                                      ##
##      https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_EVI_5km_Monthly       ##
##                                                                                                   ##
##      Downloaded directly using Google Earth Engine and hence no processing required               ##
##                                                                                                   ##
#######################################################################################################
EVI_2001 <- raster("Raw/EVI_2001.tif"); EVI_2002 <- raster("Raw/EVI_2002.tif")
EVI_2003 <- raster("Raw/EVI_2003.tif"); EVI_2004 <- raster("Raw/EVI_2004.tif")
EVI_2005 <- raster("Raw/EVI_2005.tif"); EVI_2006 <- raster("Raw/EVI_2006.tif")
EVI_2007 <- raster("Raw/EVI_2007.tif"); EVI_2008 <- raster("Raw/EVI_2008.tif")
EVI_2009 <- raster("Raw/EVI_2009.tif"); EVI_2010 <- raster("Raw/EVI_2010.tif")
EVI_2011 <- raster("Raw/EVI_2011.tif"); EVI_2012 <- raster("Raw/EVI_2012.tif")
EVI_2013 <- raster("Raw/EVI_2013.tif"); EVI_2014 <- raster("Raw/EVI_2014.tif")
EVI_2015 <- raster("Raw/EVI_2015.tif")
EVI_stack <- stack(EVI_2001, EVI_2002, EVI_2003, EVI_2004, EVI_2005, EVI_2006, EVI_2007, EVI_2008, 
                   EVI_2009, EVI_2010, EVI_2011, EVI_2012, EVI_2013, EVI_2014, EVI_2015)
EVI_mean <- mean(EVI_stack)
raster::writeRaster(EVI_mean, filename = "Processed/EVI_Overall_Mean.tif", format = "GTiff", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##          Landcover (Dominant Class (Band 1) and Individual Landtypes (Remaining Bands)            ##
##                                                                                                   ##
##      Landcover (Yearly Average 2001-2012 and Dominant Class for 2001-2012)                        ##
##        Rasters of the fractional international geosphere biosphere programme landcover            ##
##        Multiple bands (hence importing as a stack), with each band representing the %             ##
##        of each class assigned to a given pixel. First band is categorical and assigns pixels      ##
##        dominant class of each pixel. From the Malaria Atlas Project. 20 bands per stack.          ##
##                                                                                                   ##
##      List of Bands:                                                                               ##
##        1. Overall Class (Dominant class of each resulting pixel)                                  ##
##        2. Water   3. Evergreen Needleaf Forest   4. Evergreen Broadleaf Forest                    ##
##        5. Deciduous Needleleaf Forest   6. Deciduous Broadleaf Forest   7. Mixed Forest           ##
##        8. Closed Shrublands   9. Open Shrublands   10. Woody Savanas   11. Savannas               ##
##        12. Grasslands   13. Permanent Wetlands   14. Cropland   15. Urban                         ##
##        16. Cropland/Natural Veg Mosaic   17. Snow & Ice   18. Barren/Sparsely Populated           ##
##        19. Unclassified   20. No Data                                                             ##
##                                                                                                   ##
##    Only using layer 1, and taking the mode for the years 2001-2012                                ##
##                                                                                                   ##
##    https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_IGBP_Fractional_Landcover_5km_Annual 
##                                                                                                   ##
#######################################################################################################
Landcover_2001 <- stack("Raw/Landcover_2001_5000.tif"); Landcover_2002 <- stack("Raw/Landcover_2002_5000.tif")
Landcover_2003 <- stack("Raw/Landcover_2003_5000.tif"); Landcover_2004 <- stack("Raw/Landcover_2004_5000.tif")
Landcover_2005 <- stack("Raw/Landcover_2005_5000.tif"); Landcover_2006 <- stack("Raw/Landcover_2006_5000.tif")
Landcover_2007 <- stack("Raw/Landcover_2007_5000.tif"); Landcover_2008 <- stack("Raw/Landcover_2008_5000.tif")
Landcover_2009 <- stack("Raw/Landcover_2009_5000.tif"); Landcover_2010 <- stack("Raw/Landcover_2010_5000.tif")
Landcover_2011 <- stack("Raw/Landcover_2011_5000.tif"); Landcover_2012 <- stack("Raw/Landcover_2012_5000.tif")
Landcover_Stack <- stack(Landcover_2001[[1]], Landcover_2002[[1]], Landcover_2003[[1]], Landcover_2004[[1]], Landcover_2005[[1]], Landcover_2006[[1]],
                         Landcover_2007[[1]], Landcover_2008[[1]], Landcover_2009[[1]], Landcover_2010[[1]], Landcover_2011[[1]], Landcover_2012[[1]])
Dominant_Landcover <- calc(Landcover_Stack, fun = raster_mode)
writeRaster(Dominant_Landcover, filename = "Processed/Dominant_Landcover.tif", format = "GTiff", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##         Water Sources (Inland Water Areas and Rivers etc from Digital Chart of the World)         ##
##                                                                                                   ##
##    Continuous surface of straight-line distance to nearest water body generated from a map of     ##
##    the rivers, canals and lakes in India. Separate files for lines (rivers etc) and area          ##
##    (lakes etc) features.                                                                          ##
##                                                                                                   ## 
##    Provided in country size shapefiles.                                                           ##
##    http://www.diva-gis.org/gdata (Unprocessed, vector and line map of rivers, canals and lakes)   ## 
##                                                                                                   ##
#######################################################################################################
Water_Areas <- read_sf("Raw/Water_Areas_DIVA_GIS/IND_water_areas_dcw.shp")
Water_Areas_Spatial_Polygons <- as_Spatial(Water_Areas)

beginCluster()
Water_Areas_Spatial_Polygons_Lines <- as(Water_Areas_Spatial_Polygons, "SpatialLinesDataFrame")
Water_Areas_Spatial_Polygons_Lines@data[, 1] <- 1
water_areas_blank_raster_lines <- raster(nrow = 1000, ncol = 1000, xmn = 60, xmx = 100, ymn = 0, ymx = 40)
Water_Areas_Raster_lines <- rasterize(Water_Areas_Spatial_Polygons_Lines, water_areas_blank_raster_lines, background = 0) 
Water_Areas_Raster_lines[Water_Areas_Raster_lines > 0] <- 1

Water_Areas_Spatial_Polygons@data[, 1] <- 1
water_areas_blank_raster <- raster(nrow = 1000, ncol = 1000, xmn = 60, xmx = 100, ymn = 0, ymx = 40)
Water_Areas_Raster <- rasterize(Water_Areas_Spatial_Polygons, water_areas_blank_raster, background = 0) 
Water_Areas_Raster[Water_Areas_Raster > 0] <- 1
endCluster()

combined <- mosaic(Water_Areas_Raster_lines, Water_Areas_Raster, fun = max)
raster::writeRaster(combined, filename = "Processed/Water_Areas_DIVA_GIS_Processed.tif", format = "GTiff", overwrite = TRUE)

combined[combined == 0] <- NA
Distance_to_Water_DCW <- distance(combined, progress = "text")
raster::writeRaster(Distance_to_Water_DCW, filename = "Processed/DCW_Distance_to_Water.tif", format = "GTiff", overwrite = TRUE)

# LINES NOT DONE YET AS WAS TAKING AGES. MIGHT NOT BE NECESSARY
# beginCluster()
# Water_Lines <- read_sf("Raw/Water_Areas_DIVA_GIS/IND_water_lines_dcw.shp")
# Water_Lines_Spatial_Lines <- as_Spatial(Water_Lines)
# water_lines_blank_raster <- raster(nrow = 1000, ncol = 1000, xmn = 60, xmx = 100, ymn = 0, ymx = 40)
# Water_Lines_Raster <- rasterize(Water_Lines_Spatial_Lines, water_lines_blank_raster, background = 0) # consider also doing fun = count
# endCluster()
# plot(Water_Lines_Spatial_Lines)
# plot(Water_Lines_Raster, xlim = c(80, 85), ylim = c(25, 30))
# raster::writeRaster(Water_Lines_Raster, filename = "Processed/Water_Lines_DIVA_GIS_Processed.tif", format = "GTiff")


#######################################################################################################
##                                                                                                   ##
##        Distance to Nearest Water Bodies (from WWF's Global Lakes and Wetlands Database)           ##
##  Lakes, reservoirs, rivers and different wetland types. There are three layers to this data.      ##
##  Layer 1 comprises the 3067 largest lakes (area >= 50km2) and 654 largest reservoirs worldwide.   ##
##  Layer 2 comprises permanent open water bodies with a surface area of >= 0.1km2, excluding those  ##
##  contained in GLWD-1. Layer 3 combines these polygons with other info about types of wetlands to  ##
##  create a global raster map of wetland/water body extent at 30-second resolution.                 ##
##                                                                                                   ##
##  For Layer 3, the values mean the following things:                                               ##
##      1 Lake  2 Reservoir  3 River  4 Freshwater Marsh, Floodplain  5 Swamp Forest/Flooded Forest  ##
##      6 Coastal Wetland (inc Mangroves etc)  7 Pan, Brackish, Saline Wetland                       ##
##      8 Bog, Fen, Mire (Peatland)  9 Intermittent Wetland/Lake  10 50-100% Wetland                 ##
##      11 25-50% Wetland   12 Wetland Complex (0-25% Wetland)                                       ##
##                                                                                                   ##
##  https://www.worldwildlife.org/publications/global-lakes-and-wetlands-database-lakes-and-wetlands-grid-level-3
##  Note: Obi uses this data to calculate distance to nearest water bodies, and the other (Digital   ##
##        Chart of the World) to calculate distance to permanent rivers.                             ##
##                                                                                                   ##
#######################################################################################################
Small_Water_Bodies <- read_sf("Raw/GWLD_Level_2/glwd_2.shp")
Small_Water_Bodies <- as_Spatial(Small_Water_Bodies)
Small_Water_Bodies <- crop(Small_Water_Bodies, India_Surroundings_Extent)

Large_Water_Bodies <- read_sf("Raw/GWLD_Level_1/glwd_1.shp")
Large_Water_Bodies <- as_Spatial(Large_Water_Bodies)
Large_Water_Bodies <- crop(Large_Water_Bodies, India_Surroundings_Extent)
Large_Water_Bodies <- Large_Water_Bodies[, names(Small_Water_Bodies)]

All_Water_Bodies <- rbind(Large_Water_Bodies, Small_Water_Bodies)
All_Water_Bodies@data[, 1] <- 1
water_areas_blank_raster <- raster(nrow = 1000, ncol = 1000, xmn = 60, xmx = 100, ymn = 0, ymx = 40)
Water_Areas_Raster_Polygons <- rasterize(All_Water_Bodies, water_areas_blank_raster, background = 0) 
Water_Areas_Raster_Polygons[Water_Areas_Raster_Polygons > 0] <- 1

All_Water_Bodies_Lines <- as(All_Water_Bodies, "SpatialLinesDataFrame")
All_Water_Bodies_Lines@data[, 1] <- 1
water_areas_blank_raster <- raster(nrow = 1000, ncol = 1000, xmn = 60, xmx = 100, ymn = 0, ymx = 40)
Water_Areas_Raster_Lines <- rasterize(All_Water_Bodies_Lines, water_areas_blank_raster, background = 0) 
Water_Areas_Raster_Lines[Water_Areas_Raster_Lines > 0] <- 1

combined <- mosaic(Water_Areas_Raster_Polygons, Water_Areas_Raster_Lines, fun = max)
raster::writeRaster(combined, filename = "Processed/WWF_Layers_1_&_2_Processed.tif", format = "GTiff", overwrite = TRUE)

combined[combined == 0] <- NA
Distance_to_Water_WWF <- distance(combined, progress = "text")
plot(Distance_to_Water_WWF)
raster::writeRaster(Distance_to_Water_WWF, filename = "Processed/WWF_Distance_to_Water.tif", format = "GTiff", overwrite = TRUE)


#######################################################################################################
##                                                                                                   ##
##                           Max Cumulative Rainfall In Three Months                                 ##
##                                                                                                   ##
#######################################################################################################
number_to_do <- 830452 # 830452
dates <- date_generator(1981)

rainfall_stack <- matrix(nrow = number_to_do, ncol = length(dates))
for (i in 1:length(dates)) {
  file_name <- paste0("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/Datasets/CHIRPS_Rainfall_Data/Overall_India_CHIRPS_Rainfall/rainfall_", dates[i], ".tif")
  temp_raster <- raster(file_name)
  rainfall_stack[, i] <- as.vector(temp_raster)[1:number_to_do] # check whether as.vector takes the raster row by row or column by column
  print(i)
}

total_rainfall <- c()
for (i in 1:number_to_do) {
  total_rainfall[i] <- sum(rainfall_stack[i, ], na.rm = TRUE)
  if (i %% 10000 == 0) {
    print(i)
  }
}

cumulative_rainfall <- matrix(nrow = number_to_do, ncol = length(dates))
for (i in 1:number_to_do) {
  for (j in 1:length(dates))
    if (j <= length(dates) - 89) {
      cumulative_rainfall[i, j] <- sum(rainfall_stack[i, j:(j+89)], na.rm = TRUE)
    } else {
      indices_to_year_end <- seq(from = j, to = length(dates))
      start_year_stop <- 90 - length(indices_to_year_end)
      indices_start_year <- seq(1:start_year_stop)
      overall_indices <- c(indices_start_year, indices_to_year_end)
      cumulative_rainfall[i, j] <- sum(rainfall_stack[i, overall_indices], na.rm = TRUE)
    }
  if (i %% 1000 == 0) {
    print(i)
  }
}

prop_cumu_rainfall <- matrix(nrow = number_to_do, ncol = length(dates))
for (i in 1:number_to_do) {
  if (total_rainfall[i] != 0) {
    prop_cumu_rainfall[i, ] <- cumulative_rainfall[i, ]/total_rainfall[i]
  } else {
    prop_cumu_rainfall[i, ] <- 0
  }
  if (i %% 10000 == 0) {
    print(i)
  }
}

max_prop_cumu_rainfall <- apply(prop_cumu_rainfall, 1, max)
mat_version <- matrix(max_prop_cumu_rainfall, nrow = 931, ncol = 892, byrow = TRUE)

extent <- extent(temp_raster) 
r <- raster(mat_version)
extent(r) <- extent
plot(r)



#######################################################################################################
##                                                                                                   ##
##                           Other Malaria Atlas Project Covariates                                  ##
##                                                                                                   ##
##     In addition to the covariates collated above, a number of covariates used in previous         ##
##     Malaria Atlas Project publications were utilised (see https://map.ox.ac.uk/ for more          ## 
##     details).                                                                                     ##
##                                                                                                   ##
##     Downloaded directly from MAP servers hence no processing required                             ##
##                                                                                                   ##
#######################################################################################################
Wetlands <- raster("Processed/WWF_Global_Lakes_Wetlands_Level3.tif")
Wetlands <- crop(Wetlands, India_Surroundings_Extent)
Wetlands[Wetlands >= 1] <- 1 
Wetlands[Wetlands < 1] <- NA 
Distance_to_Wetlands <- distance(Wetlands, progress = "text")
plot(Distance_to_Wetlands)
raster::writeRaster(Distance_to_Wetlands, filename = "Processed/WWF_Distance_to_Wetlands.tif", format = "GTiff", overwrite = TRUE)

