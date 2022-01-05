#######################################################################################################
##                                                                                                   ##
##                        Loading Required Libraries, Shapefiles and Data                            ##
##                                                                                                   ##
#######################################################################################################
library(rgdal); library(sf); library(raster); library(sp); library(maps); library(zoo)
source("Functions/Spatial_Data_Processing_Functions.R")
location_metadata <- read.csv("Datasets/Systematic_Review/Geographical_Data.csv")
location_data <- location_metadata[, c(3, 7, 8, 12, 13, 14:37)] 


#######################################################################################################
##                                                                                                   ##
##   Processing Data Subset #1 - Point Locations I.e. Actual Points & Polygons Where Buffers         ##
##                               Are To Be Generated Around Single Pairs of Geocoordinates           ##
##                                                                                                   ##
##  Locations where coordinates with high resolution and certainty were assigned as points here. In  ##
##  other instances however, there was uncertainty in the exact location. This stems from:           ##
##                                                                                                   ## 
##     1) Uncertainty in exact sampling location.                                                    ##
##     2) The resolution of the coordinates was not high.                                            ##
##     3) Single pair of geocoordinates given covering multiple sampling sites                       ##
##                                                                                                   ##
##  In these instances, to reflect the uncertainty surrounding these geocoordinate points, the       ##
##  code below generates a buffer around the coordinate pair to create a circular polygon.           ##
##  It does this in the following way:                                                               ##
##                                                                                                   ##
##     1) Subset the dataframe into points associated with a particular polygon size and which ARE   ##
##        NOT associated with a particular administrative unit or set of coordinates that can be     ##
##        used to construct a polygon manually.                                                      ##
##     2) Creates a vector of UTM codes to reproject the coordinates into. Necessary as a            ##
##        circle drawn in one projection system won't be a proper circle in another.                 ##
##     3) Create a new sf object using each row of the vector of points classified as a particular   ##
##        polygon type (WIDE AREA, POLYGON SMALL or POLYGON LARGE).                                  ##
##     4) Transform that coordinate pair onto a linear projection using its UTM code.                ##
##     5) Create a buffer round that point in the UTM projection.                                    ##
##     6) Reproject that polygon back into the coordinate system I'm working in.                     ##
##                                                                                                   ##
##  Note: Same process applied to WIDE AREAS (circular buffer, 25km^2), SMALL POLYGONS               ##
##        (circular buffer, 100km^2) and LARGE POLYGONS (circular, buffer, 250km^2).                 ##
##                                                                                                   ##
##  Note: Some of the constructed polygons extend into the sea. Manual correction later on is to     ##
##        exclude NA covariate values arising from the sea.                                          ##
##                                                                                                   ##
#######################################################################################################

#### POINTS ####
points_raw <- location_data[location_data$Geo_Type == "Single" & location_data$Resolution == "Point", 1:7]
points <- st_as_sf(points_raw, coords = c("Lon_1", "Lat_1"), crs = '+proj=longlat +datum=WGS84')                  

#### WIDE AREAS ####
wide_areas_raw <- location_data[location_data$Geo_Type == "Single" & location_data$Resolution == "Wide Area", 1:7]
wide_areas <- generate_buffers(wide_areas_raw, "single")

#### SMALL POLYGONS ####
small_polygons_raw <- location_data[location_data$Geo_Type == "Single" & location_data$Resolution == "Small Polygon", 1:7]
small_polygons <- generate_buffers(small_polygons_raw, "single")

#### LARGE POLYGONS ####
large_polygons_raw <- location_data[location_data$Geo_Type == "Single" & location_data$Resolution == "Large Polygon", 1:7]
large_polygons <- generate_buffers(large_polygons_raw, "single")

# Combining into single dataframe
single_locations <- rbind(points, wide_areas, small_polygons, large_polygons)

########################################################################################################
##                                                                                                    ##
##               Processing Data Subset #2 - Polygons with Only 2 Geocoordinate Pairs                 ##
##                                                                                                    ##
##  The code below is executed on data where polygons are composed only of two pairs of coordinates.  ##
##  In practice, this arises as some studies were done with only 2 villages but data was presented    ##
##  for both villages aggregated together. Or else, there were cases where >2 villages were surveyed  ##
##  (and data presented aggregated together) but it was only possible to geolocate 2 villages.        ##
##  Cannot create polygon based on only two datapoints and so these instances are handled             ##
##  differently. Specifically, the covariate value for each individual geocoordinate pair will be     ##
##  extracted, and then the average of those two values will be taken.                                ##
##                                                                                                    ##
##     1) Remove any data entries which are composed of more than 2 pairs of geocoordinates. Are      ##
##        dealt with in "Processing Data Subset #3".                                                  ##
##     2) Subset the data to produce a dataset with only 2 pair polygons.                             ##
##     3) Create a small polygon around each of the geocoordinate pairs.                              ##
##     4) Add those buffered points to the overall dataframe. Will deal with covariate averaging      ##
##        downstream.                                                                                 ##
##                                                                                                    ##
##  Note: Number of instances of polygons with only two pairs of points. A polygon then can't be      ##
##        constructed. Construct a buffer round each point, then average the covariate value          ## 
##        over these?                                                                                 ##
##                                                                                                    ##
########################################################################################################
pair_locations_raw <- location_data[location_data$Geo_Type == "2_Singles", 1:9]
pair_locations <- generate_buffers(pair_locations_raw, "pair")

########################################################################################################
##                                                                                                    ##
##      Processing Data Subset #3 - Polygons to be Constructed Based on Geocoordinate Pair Sets       ##
##                                                                                                    ##
##  In a number of instances, the paper provides (or we were able to find) multiple geocoordinate     ##
##  pairs pertaining to multiple sampling sites, but then presents the data aggregated together.      ##
##  In these cases, we explicitly construct polygons based on the provided sets of points.            ##
##                                                                                                    ##
########################################################################################################
polygons_raw <- location_data[location_data$Geo_Type == "Polygon" & location_data$Resolution == "Constructed", ]
constructed_polygons <- construct_polygons(polygons_raw)

# Combining All 3 Types of Area to Generate the Overall Dataframe
location_coordinates <- rbind(single_locations, constructed_polygons, pair_locations)

