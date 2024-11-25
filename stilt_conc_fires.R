# Convolve footprints with FINNv2.5 emissions and WBSE emissions
library(dplyr)
library(raster)
library(ncdf4)
library(rgdal) # package for geospatial analysis
library(lubridate)
library(splitstackshape) #for cSplit function
library(stringi) #for getings parts of a string
library(hash)
# Set working directory
setwd("E:/Claire")
stilt_wd <- 'out_continuous'
getwd()

## # Get the daily average footprints to convolve with emissions inventories
daily_footprints <- list.files(path = file.path(stilt_wd, 'daily_avg_cont_footprints'), 
                               pattern = "nc", full.names = FALSE, recursive = FALSE)
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

# Read data and save PM2.5 emissions data
# Load emissions inventory, assigning the gridded product to "emissions" and
# extracting the time for each grid to "emissions_time"
WBSE_fires <- list.files(path=file.path('individual_fires'), pattern='.tif', full.names=FALSE, recursive = FALSE)
d_fires <- hash()
for (fire in WBSE_fires){
  emissions <- brick(file.path('individual_fires', fire)) # kg/day
  ## Convert emissions from kg/day to kg/m2/day
  area_3km_grid <- raster::area(emissions) # get the area of each cell in km^2
  area_3km_grid_meters <- area_3km_grid*1000000 # km^2 --> m^2
  emissions_kg_m2 <- emissions/area_3km_grid_meters #kg/m^2
  emissions_kg_m2 <- raster::resample(emissions_kg_m2, footprint_ex, method='ngb')
  dates <- as.data.frame(seq(from = as.Date("2018/1/1"), by="days", length.out = 365))
  emissions_kg_m2 <- setZ(emissions_kg_m2, dates[,1], "EmisDate")
  emissions_time <- getZ(emissions_kg_m2)
  .set(d_fires, keys=fire, values=emissions_kg_m2)
}

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
  concentrations <- data.frame((matrix(ncol = 16, nrow = length(string_dates))))
  colnames(concentrations) <- c('Time_PST', names(d_fires))
  concentrations$Time_PST <- STILT_dates
  # Calculate concentrations 
  for (i in 1:length(string_dates)){
    # Import footprints and extract timestamps
    daily_footprints_day <- grep(string_dates[i],daily_footprints_loc, value=TRUE) 
    for (f in 1:length(WBSE_fires)){
      fire = names(d_fires)[f]
      fire_emissions = d_fires[[fire]]
      fire_conc=0
      if(!is_empty(daily_footprints_day)){
          for (x in 1:length(daily_footprints_day)){
            footday_path <-file.path(stilt_wd,'daily_avg_cont_footprints', daily_footprints_day[x])
            foot <- brick(footday_path) # ppm/(umol m-2 s-1)
            nc_time <- as.Date(c(getZ(foot)), tz= 'PST', origin = '1970-01-01')
            
            # Subset emissions to match the footprint timestamps
            band <- findInterval(nc_time, emissions_time)
            fire_subset <- subset(fire_emissions, band)
            
            #Calculate the near-field PM2.5 contribution by taking the product 
            # of the footprints and the emissions fluxes
            # sum all concentrations from days of back trajectories
            fire_conc = sum(raster::values(foot * fire_subset * emi_kg_to_umol), na.rm = T)*1000 * wght_pm25 *(1/ 24.45) + fire_conc
          }
      }
      concentrations[i,f+1] <- fire_conc
    }
  }
  # Write csv for all STILT modeling dates 
  write.table(concentrations, file = file.path(stilt_wd, paste0(coordsString, "_individfires_conc_all2018days_72.csv")), 
              append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
}

