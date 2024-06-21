# Libraries
#install.packages(c('tidyverse', 'raster'))
library(tidyverse)
library(raster)
library(ncdf4)
library(pbapply)
library(splitstackshape)
library(lubridate)
#stilt_receptor_coordinates <- readRDS("E:/Claire/Claire_abl_24/stilt_receptor_coordinates.rds")
stilt_dates <- readRDS("E:/Claire/out_continuous/stilt_dates.rds")
footprints_expanded <- readRDS("E:/Claire/out_continuous/out_continuous_expanded.rds")

stilt_directory_path <- "E:/Claire/out_continuous/by-id/"
stilt_wd <- "E:/Claire/out_continuous/"
setwd(stilt_wd)

# Get all of the footprints from all receptors and dates
footprints <- list.dirs(path = file.path(stilt_directory_path), full.names = FALSE, recursive = FALSE)
footprints2 <- data.frame(footprints)
footprints3 <- cSplit(footprints2, "footprints", "_")
names(footprints3) <- c("run_time", "long", "lati", "zagl")
footprints3$all_dates <- substr(footprints3$run_time, 1, 8)
string_dates <- unique(footprints3$all_dates)
STILT_dates <- ymd(string_dates)

footprints3$location <- paste0(footprints3$long, "_", footprints3$lati)
loc_unique <- unique(footprints3$location)
length(loc_unique)

# Generate Daily Average Footprints for each day and receptor 
for (jj in 1:length(loc_unique)){
  message("% done: ", (jj -1) * 100.0 / length(loc_unique))
  
  coordsString <- loc_unique[jj]
  footprint_paths <- footprints_expanded$path[footprints_expanded$receptor == coordsString]
  length(footprint_paths)
  
  for (i in 1:length(STILT_dates)){
    # For each footprint in "footprint_date", fetch the footprint total (all timesteps in a day) into a list
    # Get 4 footprints for each day 
    daily_footprints <- grep(string_dates[i], footprint_paths, value=TRUE)
    if (!is_empty(daily_footprints)){
        footlist_day <- lapply(1:length(daily_footprints), function(x) {
        # Import footprint and extract timestamp
        footday_path <- file.path(daily_footprints[x])
        if (!is.null(footday_path) && file.exists(footday_path)==TRUE){
          footday <- brick(footday_path) # ppm/(umol m-2 s-1)
          footday
        }})
        footlist_day <- footlist_day %>% discard(is.null)
        # Get the list of unique dates in the four backward footprints for a given day
        # loop & append
        actual_dates_vector <- vector()
        for(y in 1:length(footlist_day)){
          dts <- as.Date(getZ(footlist_day[[y]]))
        
          # convert dts (Date class objects) to character with the desired format
          dts <- format(dts,format="%Y-%m-%d")
        
          actual_dates_vector <- c(actual_dates_vector,dts)
        }
        actual_dates_unique <- unique(as.Date(actual_dates_vector))
        for (z in 1:length(actual_dates_unique)){
          # Subset the footprint layers for each day of back trajectories 
          substacks <- lapply(1:length(footlist_day), function(a){
            selected <- which(as.Date(getZ(footlist_day[[a]])) == actual_dates_unique[[z]])
            if (length(selected) >0){
              sub <- subset(footlist_day[[a]], selected)
              substack <- stack(sub)
              substack
            }
          })
          # Generate average footprint for each day of back trajectories
          foot_average <- brick(sum(stack(Filter(Negate(is.null), substacks))) / length(Filter(Negate(is.null),substacks)))
          foot_average <- setZ(foot_average, actual_dates_unique[[z]], name='FootprintDate')
          writeRaster(foot_average, filename=file.path(stilt_wd, 'daily_avg_cont_footprints',paste0(string_dates[i],"_",coordsString, "_", actual_dates_unique[[z]], "_", "dailyavg.nc")), overwrite=TRUE, format="CDF", varname="Sensitivity", varunit="ppm (umol-1 m2 s)", 
                    longname="Sensitivity -- raster stack to netCDF", xname="Longitude",   yname="Latitude", zname="Date")
        }
      }  
    }
    
}

# Get all of the daily footprints for each receptor and create seasonal average footprint for fire season 2018 
daily_footprints <- list.files(path = file.path('daily_avg_cont_footprints'), full.names = FALSE, recursive = FALSE)
daily_footprints2 <- data.frame(daily_footprints)
daily_footprints3 <- cSplit(daily_footprints2, 'daily_footprints', '_')
names(daily_footprints3) <- c("date", "long", "lati", "time_interval", "ending")
daily_footprints3$location <- paste0(daily_footprints3$long, "_", daily_footprints3$lati)
loc_unique <- unique(daily_footprints3$location)

for (jj in 1:length(loc_unique)){
    all_days <- grep(loc_unique[jj], daily_footprints, value=TRUE)
    footlist_day <- lapply(1:length(all_days), function(x) {
      # Import footprint and extract timestamp
      footday_path <- file.path('daily_avg_cont_footprints', all_days[x])
      if (file.exists(footday_path) == TRUE){
        footday <- raster(footday_path) # ppm/(umol m-2 s-1)
        }
      })
    # Calculate the mean footprint for a receptor
    #foot_average <- brick(sum(stack(Filter(Negate(is.null), footlist_day))) / length(Filter(Negate(is.null),footlist_day)))
    stack <- stack(Filter(Negate(is.null), footlist_day))
    foot_average <- calc(stack, mean, na.rm=TRUE)
    writeRaster(foot_average, filename=file.path(stilt_wd, 'mean_footprints', paste0(loc_unique[jj], "_", "mean_footprint.nc")), overwrite=TRUE,  format="CDF", varname="Sensitivity", varunit="ppm (umol-1 m2 s)", 
                longname="Sensitivity -- raster stack to netCDF", xname="Longitude",   yname="Latitude")
}
