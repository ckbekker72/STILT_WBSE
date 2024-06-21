# Convolve footprints with FINNv2.5 emissions and WBSE emissions
library(dplyr)
library(raster)
library(ncdf4)
library(rgdal) # package for geospatial analysis
library(lubridate)
library(splitstackshape) #for cSplit function
library(stringi) #for getings parts of a string
# Set working directory
setwd("E:/Claire")
stilt_wd <- 'out_abl_test'
getwd()

#Part I: read data and save PM2.5 emissions data
# Load emisisons inventory, assigning the gridded product to "emissions" and
# extracting the time for each grid to "emissions_time"
FINN_emissions <- brick(file.path('finnv2.5_emissions_calPM25_2018_3kmresol.tif')) # kg/day/m2
dates <- as.data.frame(seq(from = as.Date("2018/1/1"), by="days", length.out = 365))
FINN_emissions <- setZ(FINN_emissions, dates[,1], "EmisDate")
emissions_time <- getZ(FINN_emissions)

# Load the WBSE emissions data
## original kg/day file
WBSE_emissions <-brick(file.path('allfires_WBSE.tif')) #kg/day
sum(raster::values(WBSE_emissions[[1:365]])) 
area_3km_grid <- raster::area(WBSE_emissions) # get the area of each cell in km^2
area_3km_grid_meters <- area_3km_grid*1000000 # km^2 --> m^2
WBSE_emissions_kg_m2 <- WBSE_emissions/area_3km_grid_meters #kg/m^2
#sum the layers
raster_sum_212 <-  raster::cellStats(WBSE_emissions_kg_m2[[212]], "sum", na.rm = TRUE)

## Modified file with kg/day/m2
WBSE_emissions_2 <- brick(file.path('allfires_WBSE_kg_m2.tif')) # kg/day/m2
dates <- as.data.frame(seq(from = as.Date("2018/1/1"), by="days", length.out = 365))
WBSE_emissions_2 <- setZ(WBSE_emissions_2, dates[,1], "EmisDate")
raster_sum2_212 <- raster::cellStats(WBSE_emissions_2[[212]], "sum", na.rm = TRUE)

# Get the daily average footprints to convolve with emissions inventories
daily_footprints <- list.files(path = file.path(stilt_wd, 'daily_avg_cont_footprints'), 
                               pattern = "nc", full.names = FALSE, recursive = FALSE)
#daily_footprints_aux <- grep('.aux', daily_footprints)
#daily_footprints <- daily_footprints[-1486]

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

# Get data from AQS and IMPROVE monitors for each day
daily_AQS_PM25 <- read.csv('dailyPM25_AQS_CA.csv')
daily_AQS_PM25$Date.Local <- as.Date(daily_AQS_PM25$Date.Local)
daily_improve_PM25 <- read.csv('improve2018data_wsites.csv')
daily_improve_PM25$SampleDate <- mdy(daily_improve_PM25$SampleDate)

# Get the station monitors and locations for each receptor 
selected_monitors <- read.csv('SELECTED_improve_aqs_distances.csv')
selected_monitors$location <- paste0(selected_monitors$Longitude, '_', selected_monitors$Latitude)

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
coordsString <- loc_unique[1]
daily_footprints_loc <- grep(coordsString, daily_footprints, value=TRUE)
daily_footprints_day <- grep(string_dates[1],daily_footprints_loc, value=TRUE) 
footprint_ex <- brick(file.path(stilt_wd,'daily_avg_cont_footprints', daily_footprints_day[1]))
FINN_emissions <- raster::resample(FINN_emissions, footprint_ex, method='ngb')
WBSE_emissions <- raster::resample(WBSE_emissions_kg_m2, footprint_ex, method='ngb')

# Loop through all receptor locations and all dates for that location 
for (jj in 1:length(loc_unique)){
  message("% done: ", (jj -1) * 100.0 / length(loc_unique))
  # Get the daily footprints for each receptor location
  coordsString <- loc_unique[16]
  daily_footprints_loc <- grep(coordsString, daily_footprints, value=TRUE)
  
  # Calculate concentrations 
  concentrations <- data.frame((matrix(ncol = 3, nrow = length(string_dates))))
  colnames(concentrations) <- c('Time_PST', 'PM25_FINN', 'PM25_WBSE')
  # Calculate concentrations 
  for (i in 1:length(string_dates)){
    # Import footprints and extract timestamps
    daily_footprints_day <- grep(string_dates[i],daily_footprints_loc, value=TRUE) 
    PM25_FINN=0
    PM25_WBSE=0
    if(!is_empty(daily_footprints_day)){
      for (x in 1:length(daily_footprints_day)){
        foot <- brick(file.path(stilt_wd,'daily_avg_cont_footprints', daily_footprints_day[x])) # ppm/(umol m-2 s-1)
        nc_time <- as.Date(c(getZ(foot)), tz= 'PST', origin = '1970-01-01')
      
        # Subset emissions to match the footprint timestamps
        band <- findInterval(nc_time, emissions_time)
        FINN_subset <- subset(FINN_emissions, band)
        WBSE_subset <- subset(WBSE_emissions, band)
      
        #Calculate the near-field PM2.5 contribution by taking the product 
        # of the footprints and the emissions fluxes
        # sum all concentrations from days of back trajectories
        PM25_FINN = sum(raster::values(foot * FINN_subset * emi_kg_to_umol), na.rm = T)*1000 * wght_pm25 *(1/ 24.45) + PM25_FINN
        PM25_WBSE = sum(raster::values(foot * WBSE_subset * emi_kg_to_umol), na.rm=T)* 1000 * wght_pm25 *(1/ 24.45) + PM25_WBSE
    }}
    concentrations[i,] <- c(STILT_dates[i], PM25_FINN, PM25_WBSE)  
  }
  # Save the concentrations for all days 
  # Write modeled concentration data for each receptor to a file
  if (jj == 1) {
    write.table(concentrations, file = file.path(stilt_wd, paste0(coordsString, "_wildfire_conc_all2018days_72.csv")), 
                append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    write.table(concentrations, file = file.path(stilt_wd, paste0(coordsString, "_wildfire_conc_all2018days_72.csv")), 
                append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)    
  }
  # Get the monitor location from the receptor list that matches coordsString
  concentrations$Time_PST <- as.POSIXlt.Date(concentrations$Time_PST)
  concentrations$location <- coordsString
  selected_monitors_coords <- selected_monitors[as.character(round(selected_monitors$Longitude,5)) == str_split(coordsString,"_")[[1]][1], ]
  
  # Subset monitor data to match the receptor location, STILT dates, and get the PM2.5 daily values
  if (selected_monitors_coords$SiteName %in% daily_AQS_PM25$Local.Site.Name == TRUE){
      aqs_monitor_data <- daily_AQS_PM25[daily_AQS_PM25$Date.Local %in% STILT_dates & daily_AQS_PM25$Local.Site.Name== selected_monitors_coords$SiteName & 
                                 daily_AQS_PM25$POC==1,]
      aqs_dates_values <- aqs_monitor_data[c('Date.Local', 'Arithmetic.Mean')]
      concentration_new <- data.frame(merge(concentrations, aqs_dates_values, by.x ='Time_PST', by.y='Date.Local'))
  }
  else if (selected_monitors_coords$SiteName %in% daily_improve_PM25$SiteName == TRUE){
    improve_monitor_data <- daily_improve_PM25[daily_improve_PM25$SampleDate %in% STILT_dates & daily_improve_PM25$SiteName == selected_monitors_coords$SiteName, ]
    improve_dates_values <- na.omit(improve_monitor_data[c('SampleDate', 'PM2.5')])
    concentration_new <- data.frame(merge(concentrations, improve_dates_values, by.x='Time_PST', by.y='SampleDate'))
  }
  
  # Write modeled and observed concentration data for each receptor to a file
  if (jj == 1) {
    write.table(concentration_new, file = file.path(stilt_wd, paste0(coordsString, "_wildfire_conc_2018days_72.csv")), 
                append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    write.table(concentration_new, file = file.path(stilt_wd, paste0(coordsString, "_wildfire_conc_2018days_72.csv")), 
                append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)    
  }
}

