#####################################################################
## Check for row duplicates in all LiDAR csv files for all UTM zones
#####################################################################
##----------------------
## TO DO
##----------------------

#### STILL TO DO ####
# - test foreach without doParallel things before
# - run before prog0 to check all files loaded
# - run similar script after prog1 with all produced files

#### SOLVED ####

#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

rm(list=ls()) # clear all variables

 
##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

nr.clusters = length(paramsGL$zones)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end

##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "spdep",
                      "spatstat",
                      "rgeos",
                      "maptools", 
                      "plyr",
                      "ggplot2",
                      "data.table",
                      "lubridate", 
                      "doParallel", 
                      "foreach"
                      # "profvis",
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

tic <- proc.time() # start clocking global time

cl <- makeCluster(nr.clusters)
registerDoParallel(cl)
## assigns to full.df the row-wise binding of the last unassigned object (dataframe) of the loop (in this case the one formed with "merge(pts9.mean.me....")
full.df <- foreach (z = 1:length(paramsGL$zones), .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
# for (z in 1:length(paramsGL$zones)) {

  zone = paramsGL$zones[z]
  print(paste('Prog0a, Check for duplicates on' ,zone))   # converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir, zone, fsep = .Platform$file.sep)

  ## load lidar veg metrics database
  infile = file.path(LOP_dir,'LOP_attr','plot_cloud_metrics_first_returns_gt2m',paste(zone,"_plot_cloud_metrics_first_returns_gt2m.csv",sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.1streturns <- read.csv(infile,head=TRUE,sep=",") ## metrics from first returns above 2m
  ## check if nr of rows is equal to nr of unique IDs
  if ( nrow(lidar.metrics.1streturns) == length(unique(lidar.metrics.1streturns$unique_id)) ) {
    counter.1streturns <- 1
  } else {
    counter.1streturns <- 0
  }

  ## load lidar biomass and inventory attributes database
  infile = file.path(LOP_dir,'LOP_attr','plot_inventory_attributes2',paste(zone,"_plot_inventory_attributes2.csv",sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.forest.attributes <- read.csv(infile,head=TRUE,sep=",") ## metrics from first returns above 2m
  ## check if nr of rows is equal to nr of unique IDs
  if ( nrow(lidar.metrics.forest.attributes) == length(unique(lidar.metrics.forest.attributes$unique_id)) ) {
    counter.forest.attributes <- 1
  } else {
    counter.forest.attributes <- 0
  }
  
  ## load lidar plot points shapefiles
  pnt_dir = file.path(LOP_dir,'LOP_attributed', fsep = .Platform$file.sep) 
  lidar <-  readOGR(dsn = pnt_dir, layer = paste(zone,"_points", sep=''))
  ## check if nr of rows is equal to nr of unique IDs
  if ( nrow(lidar@data) == length(unique(lidar@data$Unique_ID)) ) {
    counter.ID.point.shp <- 1
  } else {
    counter.ID.point.shp <- 0
  }
  ## check if there are points with same coordinates
  duplicIndic <- duplicated(lidar@coords, incomparables = FALSE)
  if ( is.na(table(duplicIndic)["TRUE"]) ) {
    counter.point.shp <- 1
  } else {
    counter.point.shp <- 0
  }
  unique.IDs.df <- as.data.frame(lidar@data$Unique_ID)
  duplic.points.ID <- as.data.frame(unique.IDs.df[duplicIndic, ])
  
  c(counter.1streturns, nrow(lidar.metrics.1streturns), length(unique(lidar.metrics.1streturns$unique_id)), 
    counter.forest.attributes, nrow(lidar.metrics.forest.attributes), length(unique(lidar.metrics.forest.attributes$unique_id)), 
    counter.ID.point.shp, nrow(lidar@data), length(unique(lidar@data$Unique_ID)), 
    counter.point.shp, nrow(duplic.points.ID))

}

stopCluster(cl)

full.df = as.data.frame(full.df)
rownames(full.df) <- paramsGL$zones
colnames(full.df) <- c("flag_1stReturns_IDs", "1stReturns_nrRows", "1stReturns_nrUniqueIDs", 
                       "flag_forestAttr_IDs", "forestAttr_nrRows", "forestAttr_nrUniqueIDs", 
                       "flag_pointsShp_IDs", "pointsShp_nrRows", "pointsShp_nrUniqueIDs", 
                       "flag_pointsShp_Coords", "pointsShp_nrDuplicateCoords")

write.csv(full.df, file = file.path(base_wkg_dir, "check_unique_IDs_coords.csv", sep = ''))

if ( sum(full.df[, "flag_1stReturns_IDs"]) == nrow(full.df) && sum(full.df[, "flag_forestAttr_IDs"]) == nrow(full.df) && sum(full.df[, "flag_pointsShp_IDs"]) == nrow(full.df) && sum(full.df[, "flag_pointsShp_Coords"]) == nrow(full.df)) {
  print("All OK!")
} else {
  print("Something's wrong!")
}

# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog0a, total elapsed time:",seconds_to_period(toc[3])))










## just in case...

#   ## load shapefile for lidar points in forest class as determined by sk_mask
#   ## load shapfiles of sample polygons
#   ## load shapefile of sample polygons center points
#   ## load shapefile for sampled lidar points in forest class with homogenous disturbed or undisturbed condition and 250m minimum spacing
#   ## this shapefile has unique_id FCID and POLY250ID dataframe columnes
#   pnt_dir = file.path(LOP_dir,'LOP_attributed', fsep = .Platform$file.sep) 
#   lidar.pts.fc <-  readOGR(dsn = pnt_dir, layer =paste(zone, "_points",sep='')) ## same layer loaded as "lidar" in prog0_lidar_plots_sampling_auto_GM.R   ## 1106975 lidar points with fc 
#   
#   lidar.sample <- readOGR(dsn = wkg_dir, layer =paste(zone, "_lidar_sample_250m", sep='')) ## 50172 lidar points within sample 5574 sample polygons
#   lidar.sample.polygons <- readOGR(dsn = wkg_dir, layer = paste(zone,"_sample_250m_poly", sep='')) # smallest polyg subset, the one with only 9 plots/polyg ## 50172 lidar points within sample 5574 sample polygons
#   #----------------  MOST IMPORTANT FILE: CENTERPOINTS FILE   ------------------------
#   lidar.sample.polygons.centerpoints <- readOGR(dsn = wkg_dir, layer = paste(zone,"_sample_250m_poly_centerpt", sep='')) ## 5574 center points USE THESE INSTEAD OF lidar.hex.centerpoints, USED TO PRODUCE SHP "UTMXXX_cpt_poly_250m_training_validation"
#   #----------------------------------------------------------------------------------
