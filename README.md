# Contribution of large wildfire events and burn severity classes to air pollution in California in 2018
This repository contains code, data, and graphics to accompany {ADD CITATION & DOI HERE ONCE PUBLISHED}. 
## Requirements
Analyses require use of Python and R. 
## Code Access
To clone this repository to your workspace: git clone https://github.com/ckbekker72/STILT_WBSE 
## Contents and Description
### 1) Receptor Selection
01_Monitor_receptors_new.ipynb: This file iterates through a selection process to choose the 16 receptors for STILT atmospheric modeling
### 2) Emissions Processing
02.1_WBSE_Emissions_Raster_Combined_500.ipynb: Processes the individual fire shapefiles for WBSE emissions to generate a multiband (365-day) raster emissions file for all fires. 
02.2_WBSE_BurnSeverity_Raster_Combined.ipynb: Processes the WBSE emission shapefiles to classify by burn severity. Generates emissions raster files for each burn severity class. 
02.3_Individual_fires_WBSE.ipynb: Generates emissions rasters for individual fires based on WBSE shapefiles. 
### 3) STILT Modeling
STILT folder contains code to conduct STILT atmospheric modeling for 16 receptors from 7/1 to 9/30/2018
Meteorological data can be accessed via https://www.ready.noaa.gov/data/archives/hrrr.v1/ and should be uploaded to the cal_hrrr data folder. 
Run the following code in the command line to conduct STILT modeling: 
Rscript Code/03_STILT/r/run_stilt_cont.R 
Repeat for run_stilt_cont2.R, run_stilt_cont3.R, and run_stilt_cont4.R. 
Alternatively, skip to the next section which also contains the processed daily average footprints.  
### 4) STILT Footprints and Concentrations
04.1_stilt_nc_to_tif.R
04.2_daily_avg_footprints_newr.R: Converts the STILT atmospheric footprints for 6-hr intervals to daily average footprints for each day (7/1-9/30) and receptor. 
04.3.1_stilt_footprints_to_conc2.R: Convolves the daily average footprints with daily FINNv2.5 and WBSE emissions to generate daily fire-derived PM2.5 concentrations. 
04.3.2_stilt_conc_fires.R: Convolves daily average footprints with individual fire emissions from WBSE. 
04.3.3_stilt_WBSE_severity_conc.R: Convolves daily average footprints with burn severity class emissions from WBSE. 
### 5) Figure Generation
05_WBSE_paper_figures 3.ipynb


