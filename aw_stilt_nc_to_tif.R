# Libraries
library(tidyverse)
library(raster)
library(ncdf4)
library(stars)
library(pbapply)

stilt_receptor_coordinates <- readRDS("E:/Aron/stilt_receptor_coordinates.rds")
stilt_dates <- readRDS("E:/Aron/stilt_dates.rds")
footprints_expanded <- readRDS("E:/Aron/footprints_expanded")

print(paste("Expected number of STILT footprints:",
            length(stilt_receptor_coordinates)*length(stilt_dates)*4))

starting_receptor = 1
n_receptors_to_attempt = length(stilt_receptor_coordinates)
starting_date = 1
n_dates_to_attempt = length(stilt_dates)
print(n_dates_to_attempt)

stilt_directory_path <- "E:/Claire/out_continuous/by-id/"

## Get stilt directories and assemble their ids
#
stilt_directories <- list.dirs(path = stilt_directory_path, full.names = TRUE, recursive = FALSE)
# 
print(paste("Actual number of STILT footprints:",length(stilt_directories)))
# 
get_ids_from_stilt_folder_names <- function(file_name){
  return(file_name %>%
            str_replace(stilt_directory_path,"") %>%
            str_replace("/",""))
   }
 
stilt_ids <- lapply(stilt_directories, get_ids_from_stilt_folder_names)

## Get list of dates from ids

get_date_from_stilt_id <- function(id){
   return(str_split(id,"_")[[1]][1] %>% 
            substr(1,8))
 }
# 
 stilt_dates <- lapply(stilt_ids, get_date_from_stilt_id) %>%
                 unlist() %>%
                 unique()
 
saveRDS(stilt_dates,"E:/Claire/out_continuous/stilt_dates.rds")

## Build a dataframe of file names with date and receptor columns
#
build_file_path <- function(stilt_id){
   return(paste0(stilt_directory_path,stilt_id,"/",stilt_id,"_foot.nc"))
 }
# 
file_paths <- lapply(stilt_ids,build_file_path)
 
footprints <- data.frame(id = unlist(stilt_ids),path = unlist(file_paths))
# 
footprints_expanded <- footprints %>% 
                     separate(id, c("date","lon","lat","junk"),"_") %>%
                     mutate(receptor = paste0(lon,"_",lat)) %>%
                     dplyr::select(receptor,date,path)
 
saveRDS(footprints_expanded,"E:/Claire/out_continuous/out_continuous_expanded.rds")

# For a specific receptor:
make_stack_for_receptor <- function(receptor_number){
  print(receptor_number)
  receptor_cords <- stilt_receptor_coordinates[[receptor_number]]
  
  footprints_for_receptor <- footprints_expanded %>%
    dplyr::filter(receptor == receptor_cords)
  
  find_average_for_date <- function(date_number){
    
    date_to_assemble <- stilt_dates[[date_number]]
    print(date_to_assemble)
    
    datetimes <- c(paste0(date_to_assemble,"0000"),
                   paste0(date_to_assemble,"0600"),
                   paste0(date_to_assemble,"1200"),
                   paste0(date_to_assemble,"1800"))
    
    footprints_for_date <- footprints_for_receptor %>%
      dplyr::filter(date %in% datetimes)
    
    full_file_names_to_join <- footprints_for_date$path
    
    # b/c 124 STILT missing, if all 4 not present, use first 4 from year
    # This is a temporary solution
    if (length(full_file_names_to_join)<4){
      print(paste("MISSING!",date_to_assemble,receptor_cords))
      full_file_names_to_join <- footprints_for_receptor$path[1:4]
    }

      # Script from CB
      
      files_to_join <- lapply(1:length(full_file_names_to_join), function(i) {
        # Import footprint and extract timestamp
        brick <- brick(full_file_names_to_join[i]) # ppm/(umol m-2 s-1)
        
        # Convert 3d brick to 2d raster if only a single timestep contains influence
        if (nlayers(brick) == 1) {
          brick <- raster(brick, layer = 1)
        } else {
          brick <- sum(brick)
        }
        return(brick)
      })
      
      return(sum(stack(files_to_join)) / length(files_to_join))
    
  }
  
  dates <- seq(starting_date,starting_date+n_dates_to_attempt-1,1)
  
  daily_rasters <- sapply(dates,find_average_for_date)

  output_file <- stack(daily_rasters)
  
  output_file_name <- paste0("E:/Aron/output/",
                             sprintf("z_stilt_%03d", receptor_number),
                             ".tif")
  
  as_stars <- st_as_stars(output_file)
  write_stars(as_stars, output_file_name)
}

receptors <- seq(starting_receptor,starting_receptor+n_receptors_to_attempt-1,1)

pbsapply(receptors, make_stack_for_receptor)