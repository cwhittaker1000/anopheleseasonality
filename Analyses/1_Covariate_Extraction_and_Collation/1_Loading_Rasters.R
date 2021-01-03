#######################################################################################################
###                                                                                                 ###
###                   Loading Required Libraries and Setting Working Directory                      ###
###                                                                                                 ###
#######################################################################################################
library(rasterVis); library(raster); library(sp); library(classInt); library(colorspace); 
library(maps); library(sp); library(rgeos); library(rgdal); library(maptools); library(dplyr);
library(tidyr); library(tmap); library(maps); library(scales); library(mapproj); library(gdalUtils); 
library(latticeExtra); library(geosphere); library(sf); library(viridis); library(gdistance); 
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/Datasets/Satellite_Covariates/Processed/")

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
Annual_Mean_Temperature <- raster("Annual_Mean_Temperature.tif")
Mean_Diurnal_Range <- raster("Mean_Diurnal_Range.tif")
Isothermality <- raster("Isothermality.tif")
Temperature_Seasonality <- raster("Temperature_Seasonality.tif")
Max_Temp_Warmest_Month <- raster("Max_Temp_Warmest_Month.tif")
Min_Temp_Coldest_Month <- raster("Min_Temp_Coldest_Month.tif")
Temp_Annual_Range <- raster("Temp_Annual_Range.tif")
Mean_Temp_Wettest_Quarter <- raster("Mean_Temp_Wettest_Quarter.tif")
Mean_Temp_Driest_Quartest <- raster("Mean_Temp_Driest_Quartest.tif")
Mean_Temp_Warmest_Quarter <- raster("Mean_Temp_Warmest_Quarter.tif")
Mean_Temp_Coldest_Quarter <- raster("Mean_Temp_Coldest_Quarter.tif")
Annual_Rain <- raster("Annual_Rain.tif")
Rain_Wettest_Month <- raster("Rain_Wettest_Month.tif")
Rain_Driest_Month <- raster("Rain_Driest_Month.tif")
Rain_Seasonality <- raster("Rain_Seasonality.tif")
Rain_Wettest_Quarter <- raster("Rain_Wettest_Quarter.tif")
Rain_Driest_Quarter <- raster("Rain_Driest_Quarter.tif")
Rain_Warmest_Quarter <- raster("Rain_Warmest_Quarter.tif")
Rain_Coldest_Quarter <- raster("Rain_Coldest_Quarter.tif")


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
PET_Yearly_Average <- raster("PET_Annual.tif")
Aridity_Yearly_Average <- raster("GAI_Annual.tif")


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
India_Pop_Density_2010 <- raster("Pop_Density_Aggregated_Tenfold.tif")


#######################################################################################################
##                                                                                                   ##
##                    Land Surface Temperature (Overall Average and SD 2001-2015)                    ##
##                                                                                                   ##
##   Land Surface Temperature (Day & Night, Overall Average and SD 2001-2015)                        ##
##      Raster maps detailing the daily and nightly land surface temperature per pixel from Google   ##
##      Earth Engine's website. Originally derived by the Malaria Atlas Project and added to GEE     ##
##      as a dataset that can be manipulated. Further details on the website. Other manipulations    ##
##      of the data possible, but presented here are overall average 2001-2015 and the standard      ##
##      deviation of the average during that period.                                                 ##
##   https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_LST_Night_5km_Monthly    ##
##   https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_LST_Day_5km_Monthly      ##
##                                                                                                   ##
##      Downloaded directly using Google Earth Engine and hence no processing required               ##
##                                                                                                   ##
#######################################################################################################
Day_LST_Mean <- raster("Day_Land_Surface_Temperature_Mean.tif")
Day_LST_SD <- raster("Day_Land_Surface_Temperature_SD.tif")
Night_LST_Mean <- raster("Night_Land_Surface_Temperature_Mean.tif")
Night_LST_SD <- raster("Night_Land_Surface_Temperature_SD.tif")


#######################################################################################################
##                                                                                                   ##
##      Tasselled Cap Wetness & Tasselled Cap Brightness (Overall Average and SD 2001-2012)          ##
##                                                                                                   ##
##    Tasselled Cap Wetness (Overall Average and SD 2001-2012)                                       ##
##        From the Malaria Atlas Project. Tasseled Cap Wetness presented as the annual mean          ##
##        from 2001-2012 and the associated Standard Deviation.                                      ##
##    https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_TCW_5km_Monthly         ##
##                                                                                                   ##
##    Tasselled Cap Brightness (Overall Average and SD 2001-2012)                                    ##
##        From the Malaria Atlas Project. Tasseled Cap Brightness presented as the annual mean       ##
##        from 2001-2012 and the associated Standard Deviation.                                      ##
##    https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_TCB_5km_Monthly         ##
##                                                                                                   ##
##    Downloaded directly using Google Earth Engine and hence no processing required                 ##
##                                                                                                   ##
#######################################################################################################
Tasseled_Cap_Wetness_Mean <- raster("Tasseled_Cap_Wetness_Mean_5000.tif")
Tasseled_Cap_Wetness_SD <- raster("Tasseled_Cap_Wetness_Sd_5000.tif")
Tasseled_Cap_Brightness_Mean <- raster("Tasseled_Cap_Brightness_Mean_5000.tif")
Tasseled_Cap_Brightness_SD <- raster("Tasseled_Cap_Brightness_Sd_5000.tif")


#######################################################################################################
##                                                                                                   ##
##                                              Elevation                                            ##
##       Elevation (90m Resolution, Assumed to Not Change Over Time)                                 ##
##          SRTM digital elevation data produced by NASA, with a resolution of 90m at the equator.   ##
##          https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003             ##
##                                                                                                   ##
##       Downloaded directly using Google Earth Engine and hence no processing required              ##
##                                                                                                   ##
#######################################################################################################
Elevation <- raster("Elevation.tif")


#######################################################################################################
##                                                                                                   ##
##                                         Humidity                                                  ##
##                                                                                                   ##
##    Specific Humidity (grams of vapour per kilogram of air)                                        ##
##        Specific humidity does not vary as the temperature or pressure of a body of air changes,   ##
##        as long as moisture isn't added or taken away. Data comes from NASA's Global Land Data     ##
##        Assimilation System (GLDAS).                                                               ##             
##    https://developers.google.com/earth-engine/datasets/catalog/NASA_GLDAS_V20_NOAH_G025_T3H       ##
##                                                                                                   ##      
##    Downloaded directly using Google Earth Engine and hence no processing required                 ##
##                                                                                                   ##
#######################################################################################################
Specific_Humidity_Mean <- raster("Specific_Humidity_mean.tif")
Specific_Humidity_SD <- raster("Specific_Humidity_sd.tif")


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
EVI_Mean <- raster("EVI_Overall_Mean.tif")


#######################################################################################################
##                                                                                                   ##
##                                      Flow Accumulation                                            ##
##                                                                                                   ##
##     Flow Accumulation (Assumed to be Largely Unchanging Over Time)                                ##
##        Part of the WWF HydroSHEDS mapping product that provides a variety of relevant             ##
##        hydrogaphic information. Flow accumulation dfines the amount of upstream area              ##
##        draining into each cell. Values range from 1 (at river sources) through to very            ##
##        large numbers at the mouths of large rivers.                                               ##
##     https://developers.google.com/earth-engine/datasets/catalog/WWF_HydroSHEDS_15ACC              ##
##                                                                                                   ##
##    Downloaded directly using Google Earth Engine and hence no processing required                 ##
##                                                                                                   ##
#######################################################################################################
Flow_Accumulation <- raster("Flow_Accumulation.tif")


#######################################################################################################
##                                                                                                   ##
##                          Water Sources (Inland Water Areas and Rivers)                            ##
##                                                                                                   ##
##      Google Earth Engine raster layer which contains information on occurrence (the frequency     ##
##      with which water was present over the 1984-2015 period), seasonality (average number of      ##
##      months of the year water is present), and max_extent (binary output with  1 if water was     ##
##      detected in a place at any point over the 1984- 2015 period). Additionally also have         ##
##      recurrence data which specifies the frequency with which water returns from year to year.    ##
##                                                                                                   ##
##      Published recently here: https://www.nature.com/articles/nature20584                         ##
##                                                                                                   ##
##      https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_0_GlobalSurfaceWater    ##
##                                                                                                   ##                                                                                                   ##
##  Note: Doesn't seem to be particularly good at picking up rivers and other thin bodies of water.  ##
##        Wonder whether this has to do with the level of spatial aggregation we're operating at.    ##
##        Or maybe my idea of how many months of the year these bodies should be available is        ##
##        warped. I don't know.                                                                      ## 
##                                                                                                   ##
#######################################################################################################
Water_Areas_Max_Extent <- raster("Water_max_extent.tif")  # binary image contianing 1 anywhere water was ever detected
Water_Areas_Seasonality <- raster("Water_seasonality.tif") # number of months where water is present
Water_Areas_Occurrence <- raster("Water_occurrence.tif")  # the frequency with which water was present
Water_Areas_Recurrence <- raster("Water_Recurrence.tif") # frequency with which water returns year to year

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
#Water_Areas_WWF <- raster("WWF_Layers_1_&_2_Processed.tif")
WWF_Distance_to_Water <- raster("WWF_Distance_to_Water.tif")
# WWF_Distance_to_Wetlands <- raster("WWF_Distance_to_Wetlands.tif")
#plot(log(WWF_Distance_to_Water))

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
#Water_Areas_DCW <- raster("Water_Areas_DIVA_GIS_Processed.tif")
DCW_Distance_to_Water <- raster("DCW_Distance_to_Water.tif")


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
Dominant_Landcover <- raster("Dominant_Landcover.tif")


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
City_Accessibility <- raster("accessibility_to_cities_2015_v1.0.tif")
CHIRPS_Max <- raster("chirps-v2-0.Synoptic.Overall.max.5km.tif")
CHIRPS_Min <- raster("chirps-v2-0.Synoptic.Overall.min.5km.tif")
CHIRPS_Mean <- raster("chirps-v2-0.Synoptic.Overall.mean.5km.tif")
WC_A0 <- raster("WorldClim_TFA_Precip_a0.Synoptic.Overall.Data.5km.mean.tif")
WC_A1 <- raster("WorldClim_TFA_Precip_a1.Synoptic.Overall.Data.5km.mean.tif")
WC_A2 <- raster("WorldClim_TFA_Precip_a2.Synoptic.Overall.Data.5km.mean.tif")
WC_A3 <- raster("WorldClim_TFA_Precip_a3.Synoptic.Overall.Data.5km.mean.tif")
WC_P0 <- raster("WorldClim_TFA_Precip_p0.Synoptic.Overall.Data.5km.mean.tif")
WC_P1 <- raster("WorldClim_TFA_Precip_p1.Synoptic.Overall.Data.5km.mean.tif")
WC_P2 <- raster("WorldClim_TFA_Precip_p2.Synoptic.Overall.Data.5km.mean.tif")
WC_P3 <- raster("WorldClim_TFA_Precip_p3.Synoptic.Overall.Data.5km.mean.tif")
Urban_Footprint <- raster("Global_Urban_Footprint_5km_PropUrban.tif")
Irrigated_Areas <- raster("Irrigated_Areas_Global_5k.tif")

# GlobCover <- raster("Globcover_v2.2_2004-2006_MAP-Mastergrids_5k.tif") # near identical to Landcover so removed
# Globcover_Wet_Areas <- raster("Wet_Areas_In_GLWD_Or_Globcover.tif") # near identical to Wetlands so removed 
