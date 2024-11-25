# Convolve footprints with WBSE burn severity emissions
library(dplyr)
library(raster)
library(ncdf4)
library(lubridate)
library(purrr)
library(splitstackshape) #for cSplit function
library(stringi) #for getings parts of a string
# Set working directory
setwd("E:/Claire")
stilt_wd <- 'out_continuous'
getwd()

## # Get the daily average footprints to convolve with emissions inventories
daily_footprints <- list.files(path = file.path(stilt_wd, 'daily_avg_cont_footprints'), 
                               pattern = ".nc", full.names = FALSE, recursive = FALSE)
# Create data frame with footprint naming data
footprints2 <- data.frame(daily_footprints)
footprints3 <- cSplit(footprints2, "daily_footprints", "_")
names(footprints3) <- c("run_time", "long", "lati", "footprint_day", "interval")
footprints3$all_dates <- substr(footprints3$run_time, 1, 8)
string_dates <- unique(footprints3$all_dates)
STILT_dates <- ymd(string_dates)

footprints3$location <- paste0(footprints3$long, "_", footprints3$lati)
loc_unique <- unique(footprints3$location)
length(loc_unique)

# Get the example footprint used to resample emissions data
daily_footprints_loc <- grep(loc_unique[1], daily_footprints, value=TRUE)
daily_footprints_day <- grep(string_dates[1],daily_footprints_loc, value=TRUE) 
footprint_ex <- brick(file.path(stilt_wd,'daily_avg_cont_footprints', daily_footprints_day[1]))

## Get the WBSE emissions by burn severity class 
WBSE_high_severity <- brick(file.path('output_burn-selected/high_severity_WBSE.tif'))
WBSE_high_severity_emissions <- WBSE_high_severity[[1:365]]
raster_highseverity_212 <-  raster::cellStats(WBSE_high_severity_emissions[[212]], "sum", na.rm = TRUE)

WBSE_moderate_severity <- brick(file.path('output_burn-selected/moderate_severity_WBSE.tif'))
WBSE_moderate_severity_emissions <- WBSE_moderate_severity[[1:365]]
raster_moderateseverity_212 <-  raster::cellStats(WBSE_moderate_severity_emissions[[212]], "sum", na.rm = TRUE)
  
WBSE_low_severity <- brick(file.path('output_burn-selected/low_severity_WBSE.tif'))
WBSE_low_severity_emissions <- WBSE_low_severity[[1:365]]
raster_lowseverity_212 <- raster::cellStats(WBSE_low_severity_emissions[[212]], "sum", na.rm = TRUE)

WBSE_grass_burn <- brick(file.path('output_burn-selected/grass_burn_WBSE.tif'))
WBSE_grass_burn_emissions <- WBSE_grass_burn[[1:365]]
raster_grassburn_212 <- raster::cellStats(WBSE_grass_burn_emissions[[212]], "sum", na.rm = TRUE)

# Function to aggregate emissions, convert to kg/m2, resample to footprint, and save with emission dates 
aggregate_and_convert <- function(WBSE_severity_emissions){
  # Aggregate from 500m to 3km (factor of 6)  
  # a new aggregated raster, sum of the values
  WBSE_severity_aggregated <- aggregate(WBSE_severity_emissions, fact=6, fun='sum')
  area_3km_grid <- raster::area(WBSE_severity_aggregated) # get the area of each cell in km^2
  area_3km_grid_meters <- area_3km_grid*1000000 # km^2 --> m^2
  WBSE_severity_kg_m2 <- WBSE_severity_aggregated/area_3km_grid_meters #kg/m^2
  # resample to grid of the footprint
  WBSE_severity_kg_m2 <- raster::resample(WBSE_severity_kg_m2, footprint_ex, method='ngb')
  # Add dates to z dimension
  dates <- as.data.frame(seq(from = as.Date("2018/1/1"), by="days", length.out = 365))
  WBSE_severity_kg_m2 <- setZ(WBSE_severity_kg_m2, dates[,1], "EmisDate")
  emissions_time <- getZ(WBSE_severity_kg_m2)
  return(WBSE_severity_kg_m2)
}
# Aggregate, Convert, and resample each emissions raster 
WBSE_high_severity_kg_m2 <- aggregate_and_convert(WBSE_high_severity_emissions)
print(raster::cellStats(WBSE_high_severity_kg_m2[[212]], "sum", na.rm = TRUE))
emissions_time <- getZ(WBSE_high_severity_kg_m2)
WBSE_moderate_severity_kg_m2 <- aggregate_and_convert(WBSE_moderate_severity_emissions)
print(raster::cellStats(WBSE_moderate_severity_kg_m2[[212]], "sum", na.rm = TRUE))
WBSE_low_severity_kg_m2 <- aggregate_and_convert(WBSE_low_severity_emissions)
print(raster::cellStats(WBSE_low_severity_kg_m2[[212]], "sum", na.rm = TRUE))
WBSE_grass_burn_kg_m2 <- aggregate_and_convert(WBSE_grass_burn_emissions)
print(raster::cellStats(WBSE_grass_burn_kg_m2[[212]], "sum", na.rm = TRUE))

# emissions kg m-2 day-1 to umol m-2 hr-1; footprint in ppm/(umol m-2 s-1) 
# 1 day = 24 hr; 1 hr = 3600 seconds; 1 kg = 1000 grams
# PM2.5 mix weight: grams/mole
# 60% carbon; 10% potassium, chlorine and calcium; 30% hydrogen, oxygen and nitrogen
wght_pm25 <- 0.6 * 12.0107 + 0.1 * ((39.0983 + 35.453 + 40.078) / 3.0) + 
  0.3 * ((1.00794 + 16.00 + 14.01)/3.0)
# A mole of a given substance is the number of grams of that substance 
# that contain 6.022 × 1023 particles (molecules) of that substance
# emissions kg m-2 day-1 to umol m-2 s-1 to match with the unit of footprint
emi_kg_to_umol <- 1000 * (1/wght_pm25) * (1/24) * (1/3600) * 1000000 #1e6 for mole to umol

# Loop through all receptor locations and all dates for that location 
for (jj in 1:length(loc_unique)){
  message("% done: ", (jj -1) * 100.0 / length(loc_unique))
  # Get the daily footprints for each receptor location
  coordsString <- loc_unique[jj]
  daily_footprints_loc <- grep(coordsString, daily_footprints, value=TRUE)
  
  # Calculate concentrations 
  concentrations <- data.frame((matrix(ncol = 5, nrow = length(string_dates))))
  colnames(concentrations) <- c('Time_PST', 'WBSE_high_severity', 'WBSE_moderate_severity', 'WBSE_low_severity', 'WBSE_grass_burn')
  concentrations$Time_PST <- STILT_dates
  # Calculate concentrations 
  for (i in 1:length(string_dates)){
    # Import footprints and extract timestamps
    daily_footprints_day <- grep(string_dates[i],daily_footprints_loc, value=TRUE) 
    WBSE_highseverity_PM=0
    WBSE_moderateseverity_PM=0
    WBSE_lowseverity_PM=0
    WBSE_grassburn_PM=0
    if(!purrr::is_empty(daily_footprints_day)){
      for (x in 1:length(daily_footprints_day)){
        foot <- brick(file.path(stilt_wd,'daily_avg_cont_footprints', daily_footprints_day[x])) # ppm/(umol m-2 s-1)
        nc_time <- as.Date(c(getZ(foot)), tz= 'PST', origin = '1970-01-01')
        
        # Subset emissions to match the footprint timestamps
        band <- findInterval(nc_time, emissions_time)
        WBSE_highseverity_subset <- subset(WBSE_high_severity_kg_m2, band)
        WBSE_moderateseverity_subset <- subset(WBSE_moderate_severity_kg_m2, band)
        WBSE_lowseverity_subset <- subset(WBSE_low_severity_kg_m2, band)
        WBSE_grass_burn_subset <- subset(WBSE_grass_burn_kg_m2, band)
        #Calculate the near-field PM2.5 contribution by taking the product 
        # of the footprints and the emissions fluxes
        # sum all concentrations from days of back trajectories
        WBSE_highseverity_PM = sum(raster::values(foot * WBSE_highseverity_subset * emi_kg_to_umol), na.rm = T)*1000 * wght_pm25 *(1/ 24.45) + WBSE_highseverity_PM
        WBSE_moderateseverity_PM = sum(raster::values(foot * WBSE_moderateseverity_subset * emi_kg_to_umol), na.rm=T)* 1000 * wght_pm25 *(1/ 24.45) + WBSE_moderateseverity_PM
        WBSE_lowseverity_PM = sum(raster::values(foot * WBSE_lowseverity_subset * emi_kg_to_umol), na.rm=T)* 1000 * wght_pm25 *(1/ 24.45) + WBSE_lowseverity_PM
        WBSE_grassburn_PM = sum(raster::values(foot * WBSE_grass_burn_subset * emi_kg_to_umol), na.rm=T)* 1000 * wght_pm25 *(1/ 24.45) + WBSE_grassburn_PM
      }}
    concentrations[i,2:5] <- c(WBSE_highseverity_PM, WBSE_moderateseverity_PM, WBSE_lowseverity_PM, WBSE_grassburn_PM)  
  }
  # Save the concentrations for all days 
  # Write modeled concentration data for each receptor to a file
  if (jj == 1) {
    write.table(concentrations, file = file.path(stilt_wd, paste0(coordsString, "_wildfireseverity_conc_all2018days_72.csv")), 
                append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    write.table(concentrations, file = file.path(stilt_wd, paste0(coordsString, "_wildfireseverity_conc_all2018days_72.csv")), 
                append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)    
  }
}
