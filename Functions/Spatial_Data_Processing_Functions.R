#######################################################################################################
###                                                                                                 ###
###                        Functions for Processing of Spatial Data                                 ###
###                                                                                                 ###
#######################################################################################################

# Function that Ascertaines the UTM Zone a Geocoordinate Pair Belong To
long2UTM <- function(lat, lon) {
  zone = (floor((lon + 180)/6) %% 60) + 1
  if (lat < 0) {
    zone = paste("WGS 84 / UTM zone ", zone, "S", sep = "")
  } else {
    zone = paste("WGS 84 / UTM zone ", zone, "N", sep = "")
  }
  return(zone)
}

# Create vector of EPSG projection codes
EPSG_projection_codes <- function(subset_location_data, single_or_pair) {
  
  EPSG <- rgdal::make_EPSG() 
  
  if (single_or_pair == "single") {
    codes <- vector(length = length(subset_location_data$Location_ID))
    for (i in 1:length(subset_location_data$Location_ID)) {
      zone <- long2UTM(subset_location_data$Lat_1[i], subset_location_data$Lon_1[i])
      codes[i] <- EPSG[grep(zone, x = EPSG$note),]$code # Search for the UTM zone code in the df of codes, EPSG, and extract the right one
    }
    return(codes)
  } 
  else if (single_or_pair == "pair") {
    counter <- 1
    codes <- c()
    for (i in 1:length(subset_location_data$Location_ID)) {
      zone_1 <- long2UTM(subset_location_data$Lat_1[i], subset_location_data$Lon_1[i])
      zone_2 <- long2UTM(subset_location_data$Lat_2[i], subset_location_data$Lon_2[i])
      codes[counter] <- EPSG[grep(zone_1, x = EPSG$note),]$code # Search for the UTM zone code in the df of codes, EPSG, and extract the right one
      counter <- counter + 1
      codes[counter] <- EPSG[grep(zone_2, x = EPSG$note),]$code # Search for the UTM zone code in the df of codes, EPSG, and extract the right one
      counter <- counter + 1
    }
    return(codes)
  }
}

# Generate Buffers 
generate_buffers <- function(subset_location_data, single_or_pair) {
  
  codes <- EPSG_projection_codes(subset_location_data, single_or_pair)
  
  if (single_or_pair == "single") {
    spatial_df <- st_as_sf(subset_location_data[1, ], coords = c("Lon_1", "Lat_1"), crs = '+proj=longlat +datum=WGS84')  # Single entry as a dummy for creation
    for (i in 1:length(subset_location_data$Location_ID)) {
      if (subset_location_data$Resolution[i] == "Wide Area") {
        buffer_size <- 2820
      } else if (subset_location_data$Resolution[i] == "Small Polygon") {
        buffer_size <- 5640
      } else if (subset_location_data$Resolution[i] == "Large Polygon") {
        buffer_size <- 8920
      }
      sf <- st_as_sf(subset_location_data[i, ], coords = c("Lon_1", "Lat_1"), crs = '+proj=longlat +datum=WGS84')
      transformed <- st_transform(sf, crs = codes[i]) # Transform coordinates onto linear projection using UTM code
      buffered <- st_buffer(transformed, dist = buffer_size) # Create buffer around coordinates
      final <- st_transform(buffered, crs = 4326) # Transforming polygon back into standard projection 
      spatial_df <- rbind(spatial_df, final) # Adding this polygon to the list of polygons
    }
    spatial_df <- spatial_df[-1, ]
    return(spatial_df)
  }
  else if (single_or_pair == "pair") {
    spatial_df <- st_as_sf(subset_location_data[1, 1:7], coords = c("Lon_1", "Lat_1"), crs = '+proj=longlat +datum=WGS84')  # Single entry as a dummy for creation
    counter <- 1
    for (i in 1:length(subset_location_data$Location_ID)) {
      if (subset_location_data$Resolution[i] == "Wide Area") {
        buffer_size <- 2820
      } else if (subset_location_data$Resolution[i] == "Small Polygon") {
        buffer_size <- 5640
      } else if (subset_location_data$Resolution[i] == "Large Polygon") {
        buffer_size <- 8920
      }
      for (j in 1:2) {
        if (j == 1) {
          sf <- st_as_sf(subset_location_data[i, 1:7], coords = c("Lon_1", "Lat_1"), crs = '+proj=longlat +datum=WGS84')
        } 
        else if (j == 2) {
          sf <- st_as_sf(subset_location_data[i, c(1:5, 8:9)], coords = c("Lon_2", "Lat_2"), crs = '+proj=longlat +datum=WGS84')
        }
        transformed <- st_transform(sf, crs = codes[counter]) # Transform coordinates onto linear projection using UTM code
        buffered <- st_buffer(transformed, dist = buffer_size) # Create buffer around coordinates
        final <- st_transform(buffered, crs = 4326) # Transforming polygon back into standard projection 
        spatial_df <- rbind(spatial_df, final) # Adding this polygon to the list of polygons
        counter <- counter + 1
      } 
    }
    spatial_df <- spatial_df[-1 , ] 
    return(spatial_df)
  }
}

# Construct Polygons
construct_polygons <- function(subset_location_data) {
  polygons_list <- list()
  crdref <- CRS('+proj=longlat +datum=WGS84')
  counter <- 1
  
  for (i in 1:length(subset_location_data$Location_ID)) {
    lat_lon_coords <- na.trim(as.numeric(subset_location_data[i, 6:29]))
    lat <- lat_lon_coords[seq(1, length(lat_lon_coords), 2)]
    lon <- lat_lon_coords[seq(2, length(lat_lon_coords), 2)]
    lat_lon <- cbind(lon, lat)
    
    points <- chull(lat_lon)
    newlatlon <- lat_lon[points, ]
    polygons_list[[counter]] <- spPolygons(newlatlon, crs = crdref)
    counter <- counter + 1
  }
  
  constructed_polygons <- st_as_sf(subset_location_data[, 1:7], coords = c("Lon_1", "Lat_1"), crs = '+proj=longlat +datum=WGS84')
  for (i in 1:length(subset_location_data$Location_ID)) {
    constructed_polygons$geometry[i] <- st_polygon(list(polygons_list[[i]]@polygons[[1]]@Polygons[[1]]@coords))
  }
  
  return(constructed_polygons)
}

get_mode <- function(v) {
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

raster_mode <- function(x) {
  ux <- unique(x)
  ux <- ux[!is.na(ux)]
  ux[which.max(tabulate(match(x, ux)))]
}

