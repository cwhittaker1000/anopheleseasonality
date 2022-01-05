# Create Vector of Dates as Strings Given a Year
date_generator <- function(year) {
  if (year < 1981) {
    year <- 1981
  }
  
  # Assessing Whether the Input is a Leap Year or Not
  leap <- year %in% c(1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2024, 
                      2028, 2032, 2036, 2040, 2044, 2048, 2052, 2056, 2060, 2064, 2068, 2072, 
                      2076, 2080, 2084, 2088, 2092, 2096)
  
  # Creating a Months Vector and then Adjuting February Depending On Whether a Leap Year or Not
  months_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  months_length[2] <- months_length[2] + leap
  
  # Creates the Character Strings Required to Query the Database for Specific Days
  years <- rep(year, sum(months_length)) # Generates the Year Part of the Date
  months <- rep(1:12, months_length) # Generates the Month Part of the Date
  months[months < 10] <- paste0("0", months[months < 10]) # Converts Single Digit Month Numbers to Double Digits 
  days <- unlist(sapply(months_length, function(x){1:x})) # Generates the Day part of the Date
  days[days < 10] <- paste0("0", days[days < 10]) # Converts the Single Digit Day Numbers to Double Digits
  dates <- paste0(years, "-", months, "-", days)
  return(dates)
}

# Creates Raster of Relevant Daily Rainfall Using Vector of Strings
raster_creator <- function(dates, directory) {
  rainfall_stack <- stack()
  for (i in 1:length(dates)) {
    daily_rainfall <- paste0(directory, "/rainfall_", dates[i], ".tif", sep = "")
    rainfall_stack <- stack(rainfall_stack, daily_rainfall)
  }
  return(rainfall_stack)
}

# Extracts Rainfall for a Specific Pair of Coordinates/Set of Coordinates
# need to check the projections because when I plot a single rainfall raster stack it's looking curved 
# and unsure whether I want a flat projection (coordinates are lat/lon in decimal)
extract_daily_rainfall <- function(raw_coordinates, rainfall_stack) {
  
  if(length(raw_coordinates) > 2) {
    temp_polygon <- spPolygons(raw_coordinates)
  } else {
    temp_polygon <- matrix(raw_coordinates, ncol = 2)
    
  }
  rainfall_time_series <- c()
  for (i in 1:nlayers(rainfall_stack)) {
    temp_rainfall_single_day <- mean(unlist(raster::extract(rainfall_stack[[i]], temp_polygon)), na.rm = TRUE)
    rainfall_time_series <- c(rainfall_time_series, temp_rainfall_single_day)
  }
  return(rainfall_time_series)
}

# Combines the Above Three Functions
rain_time_series_generator <- function(years, location_ID, location_coordinates, directory) {
  
  # Generate Vector of Dates From Which to Extract Rainfall Data 
  dates <- c()
  if (length(years) == 1) {
    dates <- date_generator(years)
  } else if (length(years) > 1) {
    for (i in 1:length(years)) {
      temp_dates <- date_generator(years[i])
      dates <- c(dates, temp_dates)
    }
  }
  
  # Creates the Stacked Raster of All Relevant Daily Rainfall Rasters
  rainfall_raster_stack <- raster_creator(dates, directory) 
  
  # Generate Coordinates
  raw_coordinates <- st_coordinates(location_coordinates$geometry[location_ID])[, 1:2]
  
  # Creates the Time Series of Rainfall At a Given Location, For the Years Specified
  rainfall_time_series <- extract_daily_rainfall(raw_coordinates, rainfall_raster_stack)
  return(list(dates = dates, rainfall_time_series = rainfall_time_series))
  
}


# Load In Specific Rainfall Dataset
rainfall_data_loader <- function(single_record_dataframe) {
  base_path <- "C:/Users/cw1716/Documents/Q_Drive_Copy/PhD/Chapter 2 - Statistical Analysis Seasonal Patterns/anopheleseasonality/Datasets/CHIRPS_Rainfall_Data/Location_Specific_CHIRPS_Rainfall/"
  Loc_ID <- single_record_dataframe$Location_ID[1]
  Ref_ID <- single_record_dataframe$Ref_ID[1]
  specific_path <- paste0("LocID_", Loc_ID, "_", "RefID_", Ref_ID, ".csv")
  overall_path <- paste0(base_path, specific_path)
  rainfall_dataset <- read.csv(file = overall_path)
  return(rainfall_dataset)
}


# Generate Dates to Subset Rainfall By (Note for Particle Filter Work, I Go Back Further In Terms of Months)
# Just want to calculate correlations here so offset not required
generate_rainfall_vector <- function(single_record_dataframe) {
  
  rainfall_dataset <- rainfall_data_loader(single_record_dataframe)
  rainfall_dataset[, 1] <- as.character(rainfall_dataset[, 1])
  
  leap_years <- c(1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016)
  months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
  months_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  starting_day <- "01"
  starting_month <- "01"
  starting_year <- single_record_dataframe$Start
  if (starting_year < 1981) {
    starting_year <- 1981
    pre_1981 <- 1
  }
  starting_year <- as.character(starting_year)
  
  ending_day <- "31"
  ending_month <- "12"
  ending_year <- single_record_dataframe$End
  if (ending_year < 1981) {
    ending_year <- 1981
    pre_1981 <- 1
  }
  ending_year <- as.character(ending_year)
  
  start_index <- paste0(starting_year, "-", starting_month, "-", starting_day)
  start_index_position <- max(which(rainfall_dataset$Date == start_index))
  
  end_index <- paste0(ending_year, "-", ending_month, "-", ending_day)
  end_index_position <- min(which(rainfall_dataset$Date == end_index))
  
  if (ending_year == 1981 & starting_year == 1981) {
    print(paste0("Location ID ", single_record_dataframe$Location_ID, " end_index_position < start_index_position"))
    end_index_position <- start_index_position + 364
  }
  
  rainfall <- rainfall_dataset[start_index_position:end_index_position, ]
  if (length(rainfall) == 365 | length(rainfall) == 366) {
    rainfall <- cbind(rainfall, stringr::str_split_fixed(rainfall$Date, "-", 2))
    colnames(rainfall) <- c("Date", "Rainfall", "Year", "Month_Day")
    rainfall <- rainfall[, c("Month_Day", "Rainfall")]
    rainfall$Month_Day <- as.character(rainfall$Month_Day)
    return(rainfall)  
  } else {
    rainfall <- cbind(rainfall, stringr::str_split_fixed(rainfall$Date, "-", 2))
    colnames(rainfall) <- c("Date", "Rainfall", "Year", "Month_Day")
    averaged_rainfall <- aggregate(Rainfall ~ Month_Day, rainfall, mean)
    averaged_rainfall$Month_Day <- as.character(averaged_rainfall$Month_Day)
    return(averaged_rainfall)
  }
}

calculate_monthly_rainfall_totals <- function(single_record_dataframe, rainfall) {
  
  leap_years <- c(1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016) # rainfall only went back to 1981 - any dates before that just used 1981
  months_length <- c(10, 10, 11, 10, 9, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
  leap_year_months_length <- c(10, 10, 11, 10, 10, 9, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 11, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 11)
  
  years_considered <- seq(single_record_dataframe$Start, single_record_dataframe$End)
  if (sum(years_considered %in% leap_years) >= 1) {
    leap_year <- TRUE
  } else {
    leap_year <- FALSE
  }
  
  monthly_mean_rainfall_totals <- c()
  monthly_sum_rainfall_totals <- c()
  
  counter <- 1
  if (leap_year) {
    for (i in 1:length(leap_year_months_length)) {
      temp_mean <- mean(rainfall[counter:(counter + leap_year_months_length[i] - 1), 2])
      monthly_mean_rainfall_totals[i] <- temp_mean
      
      temp_sum <- sum(rainfall[counter:(counter + leap_year_months_length[i] - 1), 2])
      monthly_sum_rainfall_totals[i] <- temp_sum
      
      counter <- counter + leap_year_months_length[i]
    }
  } else {
    for (i in 1:length(months_length)) {
      temp_mean <- mean(rainfall[counter:(counter + months_length[i] - 1), 2])
      monthly_mean_rainfall_totals[i] <- temp_mean
      
      temp_sum <- sum(rainfall[counter:(counter + months_length[i] - 1), 2])
      monthly_sum_rainfall_totals[i] <- temp_sum
      
      counter <- counter + months_length[i]
    }
  }
  return(list(mean_monthly = monthly_mean_rainfall_totals, sum_monthly = monthly_sum_rainfall_totals))
}