covariates <- readRDS("Datasets/Extracted_Covariates_for_Modelling/Extracted_Covariates.rds")
number_covariates <- dim(covariates)[2] - 1
covariates <- covariates[order(covariates[, "location_IDs"]), ]
covariates <- covariates[, -1]
covariates <- as.data.frame(covariates)
covariates$Dominant_Landcover <- as.factor(covariates$Dominant_Landcover)
covariates <- as.matrix(one_hot(as.data.table(covariates), cols = "Dominant_Landcover"))
number_covariates <- dim(covariates)[2] 


par(oma = c(6, 6, 6, 6))
correlation <- melt(cor(covariates[, -1]))
ggplot(correlation, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

colnames(covariates)

temp <- c("Annual_Mean_Temperature", "Mean_Diurnal_Range", "Isothermality", "Temperature_Seasonality", 
          "Max_Temp_Warmest_Month", "Min_Temp_Coldest_Month", "Temp_Annual_Range", "Mean_Temp_Wettest_Quarter", 
          "Mean_Temp_Driest_Quartest", "Mean_Temp_Warmest_Quarter", "Mean_Temp_Coldest_Quarter",
          "Day_LST_Mean", "Day_LST_SD", "Night_LST_Mean", "Night_LST_SD")
temp_cov <- covariates[, temp]
temp_melt <- melt(cor(temp_cov))
ggplot(temp_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

temp_red <- c("Temperature_Seasonality", "Annual_Mean_Temperature", "Day_LST_SD", "Mean_Temp_Driest_Quartest")
temp_red_cov <- covariates[, temp_red]
cor(temp_red_cov)

rain <- c("Annual_Rain", "Rain_Wettest_Month", "Rain_Driest_Month", "Rain_Seasonality", "Rain_Wettest_Quarter", 
          "Rain_Driest_Quarter", "Rain_Warmest_Quarter", "Rain_Coldest_Quarter", "CHIRPS_Max", "CHIRPS_Min",
          "CHIRPS_Mean", "WC_A0", "WC_A1", "WC_A2", "WC_A3", "WC_P0", "WC_P1", "WC_P2", "WC_P3")
rain_cov <- covariates[, rain]
rain_melt <- melt(cor(rain_cov))
ggplot(rain_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

rain_red <- c("Annual_Rain", "Rain_Seasonality", "WC_A3", "WC_P0", "WC_P3", "CHIRPS_Min", "Rain_Coldest_Quarter")
rain_red_cov <- covariates[, rain_red]
cor(rain_red_cov)

arid <- c("PET_Yearly_Average", "Aridity_Yearly_Average", "Specific_Humidity_Mean", "Specific_Humidity_SD",
          "Tasseled_Cap_Wetness_Mean", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_Mean", "Tasseled_Cap_Brightness_SD")
arid_cov <- covariates[, arid]
arid_melt <- melt(cor(arid_cov))
ggplot(arid_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

arid_red <- c("PET_Yearly_Average", "Specific_Humidity_SD", "Tasseled_Cap_Wetness_SD", "Tasseled_Cap_Brightness_SD")
arid_red_cov <- covariates[, arid_red]
cor(arid_red_cov)

hydro <- c("Water_Areas_Max_Extent", "Water_Areas_Seasonality", "Water_Areas_Occurrence", "Water_Areas_Recurrence", 
           "DCW_Distance_to_Water", "WWF_Distance_to_Water", "Elevation", "Flow_Accumulation", "Irrigated_Areas")
hydro_cov <- covariates[, hydro]
hydro_melt <- melt(cor(hydro_cov))
ggplot(hydro_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

hydro_red <- c("WWF_Distance_to_Water", "Water_Areas_Recurrence", "Water_Areas_Seasonality", "Irrigated_Areas", "Flow_Accumulation", "Elevation")
hydro_red_cov <- covariates[, hydro_red]
cor(hydro_red_cov)


landcover <- c("Urban_Footprint", "EVI_Mean", "India_Pop_Density_2010", "City_Accessibility", "Dominant_Landcover_0", "Dominant_Landcover_2",
               "Dominant_Landcover_5", "Dominant_Landcover_7", "Dominant_Landcover_8", "Dominant_Landcover_10", 
               "Dominant_Landcover_11", "Dominant_Landcover_12", "Dominant_Landcover_13", "Dominant_Landcover_14", "Dominant_Landcover_16")
landcover_cov <- covariates[, landcover]
landcover_melt <- melt(cor(landcover_cov))
ggplot(landcover_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab ("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

landcover_red <- c("City_Accessibility", "EVI_Mean", "Dominant_Landcover_0", "Dominant_Landcover_2",
                   "Dominant_Landcover_5", "Dominant_Landcover_7", "Dominant_Landcover_8", "Dominant_Landcover_10", 
                   "Dominant_Landcover_11", "Dominant_Landcover_12", "Dominant_Landcover_13", "Dominant_Landcover_14", "Dominant_Landcover_16")
landcover_red_cov <- covariates[, landcover_red]
cor(landcover_red_cov)


reduced_subset <- cbind(temp_red_cov, rain_red_cov, arid_red_cov, hydro_red_cov, landcover_red_cov)

cor(reduced_subset)
