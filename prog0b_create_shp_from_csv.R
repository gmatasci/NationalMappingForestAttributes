## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog0b_create_shp_from_csv.R                      
## Objective: Filter out -9999 values in the attributes, create shp from CSV files and save new CSV files (lidar and forest attributes) with new unique_id for new LiDAR data from BC

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - use foreach
# - check parameters

## SOLVED:
# -V test foreach without doParallel things before -- it works
# -V run before prog0 to check all files loaded -- done, shapefiles and csvs are now ok

#### INIT --------------------------------------------------------------------

print('Prog0b: ')

rm(list=ls()) # clear all variables

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

params0b <- list()
params0b$transect.names <- c("AFRF", "CR04", "CR08", "MK10", "Quesnel", "Tofino")
# params0b$transect.names <- c("CR04", "CR08")

params0b$targ.names.plots.BOREAL <- c("elev_mean", "elev_stddev", "elev_cv", "elev_p95", "percentage_first_returns_above_2m", "percentage_first_returns_above_mean")  ## target variable names differ between BOREAL and NONBOREAL, so we have distinct lists 
params0b$targ.names.plots.NONBOREAL <- c("elev_mean", "elev_sd", "elev_cv", "elev_p95", "cover_2m", "cover_mean")

params0b$save.attr <- F   ## default setting to save only unique_id field (attributes stored in separate CSV)
# params0b$save.attr <- T   ## setting to save shps w attributes for visual checkup


#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "spdep",
                      "spatstat",
                      "rgeos",
                      "maptools", 
                      "plyr",
                      "dplyr",
                      "ggplot2",
                      "data.table",
                      "lubridate", 
                      "doParallel", 
                      "foreach"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}


#### START ------------------------------------------------------------------

tic <- proc.time() # start clocking global time


shp.reference <-  readOGR(dsn=file.path(LOP_dir, 'LOP_transects', fsep=.Platform$file.sep), layer = "UTM10S_trsct")
CRS.reference <- crs(shp.reference)


for (t in 1:length(params0b$transect.names)) {

  transect.name <- params0b$transect.names[t]
  
#### READ & CHECK FILES -----------------------------------------------------
  
  ## read lidar veg metrics database (metrics from first returns above 2m)
  infile <- file.path(LOP_dir, 'LOP_attr_BC', 'raw_csvs', paste(transect.name, "_lidarTable.csv",sep = ''), fsep=.Platform$file.sep)
  lidar.metrics.1streturns.raw <- read.csv(infile, head=TRUE, sep=",")
  if ( !nrow(lidar.metrics.1streturns.raw) == length(unique(lidar.metrics.1streturns.raw$pix_id)) ) {
    stop(sprintf("Duplicate rows in file %s_lidarTable.csv", transect.name))
  }
  
  ## filter out rows with: -9999 in elev_p95 (all the columns except cover_2m) or 0 in cover_2m 
  lidar.metrics.1streturns <- filter(lidar.metrics.1streturns.raw, elev_p95 != -9999 & cover_2m != 0)
  
  lidar.metrics.1streturns$pix_id <- sprintf("%s_%d", transect.name, lidar.metrics.1streturns$pix_id)
  colnames(lidar.metrics.1streturns)[1] <- "unique_id"
  setnames(lidar.metrics.1streturns, params0b$targ.names.plots.NONBOREAL, params0b$targ.names.plots.BOREAL)  ## change names to BOREAL dataset names
  write.csv(lidar.metrics.1streturns, file=file.path(LOP_dir, 'LOP_attr_BC', 'plot_cloud_metrics_first_returns_gt2m', paste(transect.name, "_lidarTable_uniqueID.csv",sep = ''), fsep=.Platform$file.sep), row.names=FALSE)
  
  ## read forest attributes database (derived by linear regression from LiDAR)
  infile <- file.path(LOP_dir, 'LOP_attr_BC', 'raw_csvs', paste(transect.name, "_attributeTable.csv",sep = ''), fsep=.Platform$file.sep)
  forest.attributes.raw <- read.csv(infile, head=TRUE, sep=",")
  if ( !nrow(forest.attributes.raw) == length(unique(forest.attributes.raw$unique_id)) ) {
    stop(sprintf("Duplicate rows in file %s_lidarTable.csv", transect.name))
  }
  
  ## change unique_id field of forest.attributes.raw to match format of lidar.metrics.1streturns and filter it to keep the same plots in both tables
  forest.attributes.raw$unique_id <- sprintf("%s_%d", transect.name, forest.attributes.raw$unique_id)
  forest.attributes <- filter(forest.attributes.raw, unique_id %in% lidar.metrics.1streturns$unique_id)
  write.csv(forest.attributes, file=file.path(LOP_dir, 'LOP_attr_BC', 'plot_inventory_attributes2', paste(transect.name, "_attributeTable_uniqueID.csv",sep = ''), fsep=.Platform$file.sep), row.names=FALSE)
  
  if (params0b$save.attr) {
    lidar.points <- SpatialPointsDataFrame(lidar.metrics.1streturns[, c("utm_x", "utm_y")],
                                                  lidar.metrics.1streturns,    ## the R object to convert
                                                  proj4string = CRS.reference)
  } else {
    lidar.points <- SpatialPointsDataFrame(lidar.metrics.1streturns[, c("utm_x", "utm_y")],
                                        lidar.metrics.1streturns[, c("unique_id", "utm_x", "utm_y")],    ## the R object to convert
                                        proj4string = CRS.reference)
  }
  
  duplicIndic <- duplicated(lidar.points@coords, incomparables = FALSE)
  if ( !table(duplicIndic)["FALSE"] == length(lidar.points)) {
    stop(sprintf("Duplicate coordinates in transect %s", transect.name)) 
  }
  
  if (params0b$save.attr) {
    writeOGR(lidar.points , file.path(LOP_dir, 'LOP_attributed_BC', 'WithAttributes', fsep=.Platform$file.sep) , paste(transect.name, "_points_woSumFields", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## save it for future checks
  } else {
    writeOGR(lidar.points , file.path(LOP_dir, 'LOP_attributed_BC', 'WithoutSumFields', fsep=.Platform$file.sep) , paste(transect.name, "_points_woSumFields", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## save it for future checks
  }
    
} 


# #### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog0b, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)