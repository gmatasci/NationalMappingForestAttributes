########################################################
## Saskatchewan lidar training and validation plot selection
## Merge lidar metrics and derived forest attributes 
## Remove plots will high variablilty of vegeation height to avoid plots with mixed conditions
## Generate shapefile of 75 m sqaure polygons of training and validation plots
## Generate shapefile of center point of 75 m sqaure polygons of training and validation plots
## Generate shapefile of all 9 points in 5 m sqaure polygons of training and validation plots
## Calculate mean lidar vegetation metrics and derived forest attributes for all training and validation plots
#######################################################
##----------------------
## TO DO
##----------------------

## STILL TO DO
# - where shall I overwrite the new shapefiles from "ForGiona.zip" (use these as there was an issue wrt nr. LiDAR plots >= nr. plots with LiDAR attr. >= nr. plots with forest attr.)?
# - TRN/VAL center of polygon shapefiles match with those produced by prog0, the shapefile with all the points written at end  
#   (pts9.poly.training.validation.temp.5 as _pts9_poly_250m_training_validation) does not match!!! TO CHECK AGAIN
# - ############# = TODEL
# - check sizes of shp vs csv
# -DOUG elev_cv is ok for in-plot filtering or we need CV on 3x3 polyg ?

## SOLVED
# -V UTM11S_plot_cloud_metrics_first_returns_gt2m.csv is only 7000 KB in size and Unique_IDs do not match -- with new data from Geordie now it is OK, 7000 Kb bc is a very small transect

##----------------------
## READS  
##----------------------
# - "<UTMzone>_plot_cloud_metrics_first_returns_gt2m.csv": csv table with lidar veg metrics database
# - "<UTMzone>_plot_inventory_attributes2.csv": csv table with lidar biomass and inventory attributes database
# - "<UTMzone>_points.shp": point shp with LiDAR 25x25 m plots
# - "<UTMzone>_lidar_sample_250m.shp": (from prog0) point shp with the LiDAR plots (3x3 = 9 plots or just 1 plot) ultimately selected
# - "<UTMzone>_sample_250m_poly.shp": (from prog0) polygon shp with sampling polygons covering plots' squares, after subsetting wrt forest/dist/undist
# - "<UTMzone>_sample_250m_poly_centerpt.shp": (from prog0) point shp with center of hexagons used for sampling (coincides with "<UTMzone>_lidar_sample_250m.shp" if only central plot is kept), after subsetting wrt forest/dist/undist 

##----------------------
## WRITES
##----------------------
# - "<UTMzone>_poly_250m_training_validation.shp": (for prog2) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)
# - "<UTMzone>_cpt_poly_250m_training_validation.shp": (for prog2) point shp with centerpoints of selected polygons 
# - "<UTMzone>_pts9_poly_250m_training_validation.shp": point shp with all 9 plots within selected polygons (just to visually check 9 plots configurations, not used anymore later on)
# - "lidar_metrics_mean_training_validation.csv": csv table with all observed LiDAR and forest attributes (Y) for selected TRN and VAL samples (3x3 polygons with average plot values) for all UTM zones


#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

rm(list=ls()) # clear all variables

 
##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

# ############### TO DEL
# base_wkg_dir <- 'D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg_PROG0_FOREACH'
# paramsGL$zones <- c('UTM12S', 'UTM13S', 'UTM14S')
# ############### TO DEL


##----------------------------
## SCRIPT SPECIFIC PARAMETERS
##----------------------------
params1 <- list()
params1$elev.thresh <- 30   ## threshold on elev_p95 to eliminate birds and towers
params1$CV.thresh <- 50   ## threshold on coefficient of variation to eliminate plot with too much variability
params1$trn.pct <- 0.75   ## proportion of training samples in the trn/val split

param_file_prog1 = file.path(base_wkg_dir, 'AllUTMzones_params1.Rdata', fsep = .Platform$file.sep) 
save(params1, file = param_file_prog1)

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
  cmd=sprintf('library(%s)', pack)
  eval(parse(text=cmd))
}

tic <- proc.time() # start clocking global time

cl <- makeCluster(nr.clusters)
registerDoParallel(cl)
getDoParWorkers()
## assigns to full.df the row-wise binding of the last unassigned object (dataframe) of the loop (in this case the one formed with "merge(pts9.mean.me....")
full.df <- foreach (z = 1:length(paramsGL$zones), .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
# for (z in 1:length(paramsGL$zones)) {

  # temp.tic <- proc.time() # start clocking time for each UTM zone
  
  zone = paramsGL$zones[z]
  print(paste('Prog1, TRN/VAL splitting on' ,zone))   # converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir, zone, fsep = .Platform$file.sep)
  results_dir = file.path(base_results_dir, zone,fsep = .Platform$file.sep)

  setwd(wkg_dir)

  ## load lidar veg metrics database
  ## load lidar biomass and inventory attributes database
  infile = file.path(LOP_dir,'LOP_attr','plot_cloud_metrics_first_returns_gt2m',paste(zone,"_plot_cloud_metrics_first_returns_gt2m.csv",sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.1streturns <- read.csv(infile,head=TRUE,sep=",") ## metrics from first returns above 2m
  if ( !nrow(lidar.metrics.1streturns) == length(unique(lidar.metrics.1streturns$unique_id)) ) {
    stop(sprintf("Duplicate rows in file %s_plot_cloud_metrics_first_returns_gt2m.csv", zone))
  }

  infile = file.path(LOP_dir,'LOP_attr','plot_inventory_attributes2',paste(zone,"_plot_inventory_attributes2.csv",sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.forest.attributes <- read.csv(infile,head=TRUE,sep=",") ## metrics from first returns above 2m
  if ( !nrow(lidar.metrics.forest.attributes) == length(unique(lidar.metrics.forest.attributes$unique_id)) ) {
    stop(sprintf("Duplicate rows in file %s_plot_inventory_attributes2.csv", zone))
  }
  
  ## load shapefile for lidar points in forest class as determined by sk_mask
  ## load shapfiles of sample polygons
  ## load shapefile of sample polygons center points
  ## load shapefile for sampled lidar points in forest class with homogenous disturbed or undisturbed condition and 250m minimum spacing
  ## this shapefile has unique_id FCID and POLY250ID dataframe columnes
  pnt_dir = file.path(LOP_dir,'LOP_attributed', fsep = .Platform$file.sep) 
  lidar.pts.fc <-  readOGR(dsn = pnt_dir, layer =paste(zone, "_points",sep='')) ## same layer loaded as "lidar" in prog0_lidar_plots_sampling_auto_GM.R   ## 1106975 lidar points with fc 
  
  lidar.sample <- readOGR(dsn = wkg_dir, layer =paste(zone, "_lidar_sample_250m", sep='')) ## 50172 lidar points within sample 5574 sample polygons
  lidar.sample.polygons <- readOGR(dsn = wkg_dir, layer = paste(zone,"_sample_250m_poly", sep='')) # smallest polyg subset, the one with only 9 plots/polyg ## 50172 lidar points within sample 5574 sample polygons
  #----------------  MOST IMPORTANT FILE: CENTERPOINTS FILE   ------------------------
  lidar.sample.polygons.centerpoints <- readOGR(dsn = wkg_dir, layer = paste(zone,"_sample_250m_poly_centerpt", sep='')) ## 5574 center points USE THESE INSTEAD OF lidar.hex.centerpoints, USED TO PRODUCE SHP "UTMXXX_cpt_poly_250m_training_validation"
  #----------------------------------------------------------------------------------
  ############### TO DEL
  ############## lidar.hex.centerpoints <- readOGR(dsn = wkg_dir, layer = paste(zone,"_HexPts_250m", sep='')) # largest polyg subset, the one with all polygs, USED TO PRODUCE SHP "UTMXXX_sample_250m_poly_centerpt"
  ############### also all hex centerpoints snapped yo lidat pointds
  
  ## check names dim and length of lidar datasets
  names(lidar.metrics.1streturns)
  names(lidar.metrics.forest.attributes)
  names(lidar.pts.fc)
  names(lidar.sample)
  names(lidar.sample.polygons)
  #names(lidar.sample.polygons.centerpoints)
  
  dim(lidar.metrics.1streturns) #1219548 lidar plots 73 variables
  dim(lidar.metrics.forest.attributes) #879497 lidar plots 13 variables
  length(lidar.pts.fc) #1106975 lidar plots in forest class
  length(lidar.sample) #50172 lidar plots in forest class
  length(lidar.sample.polygons) #5574 lidar plots in forest class
  # length(lidar.sample.polygons.centerpoints) #5574 sample centerpoints
  
  ############# keep only hex center points within candidate sample polygons
  ############# point FCID and POLY250ID are the same b/c POLY250ID is based on FCID of centerpoint
  ############# remove dist and which variables
  ############# export as shapefile to visually check against lidar plots
  ############# TODEL as already computed in prog0 script
  ############# lidar.sample.polygons.centerpoint <- lidar.hex.centerpoints[!is.na(over(lidar.hex.centerpoints,as(lidar.sample.polygons,"SpatialPolygons"))),]
  ############# lidar.sample.polygons.centerpoint$POLY250ID <- lidar.sample.polygons.centerpoint$FCID
  ############# lidar.sample.polygons.centerpoint.2 <- lidar.sample.polygons.centerpoint[,(3:4)]
  
  ############# lidar.sample.polygons.centerpoints.2 <- lidar.sample.polygons.centerpoints[,(5:6)]   # should be equivalent to this

  ############# TODEL as we already have saved this info in UTM13S_sample_250m_poly_centerpt.shp, the following line is based on the "hex" centerpoints (those that are not subsampled wrt forest/dist/undist)
  ############# writeOGR(lidar.sample.polygons.centerpoint.2 , wkg_dir, paste(zone, "_sample_250m_poly_centerpt", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## large number of plots have lidar metrics but no derived forest attributes
  ## this is in part becuase only 1106975 of 1219548 plots are in forest condtion class as detemined by sk_mask
  ## however that still leaves 227478 lidar plots missing derived forest attributes that need to be removed from the lidar plot sample
  ## remove lidar plots were elev-p95 is greater than 30m indicating birds and towers
  ## merge lidar metrics and forest attributes by unique_id keeping only the plots with forest attributes
  ## then remove original two data frames b/c they take up too much space
  lidar.metrics.1streturns <- subset(lidar.metrics.1streturns, elev_p95 <= params1$elev.thresh) 
  lidar.metrics.attributes <- merge(lidar.metrics.1streturns, lidar.metrics.forest.attributes, by.x =  "unique_id", by.y = "unique_id")
  rm(lidar.metrics.1streturns)
  rm(lidar.metrics.forest.attributes)
  gc()
  dim(lidar.metrics.attributes) #879482 plots 85 variables
  

  #######################################################################
  ## next need to figure out how many of the lidar sample plots have 9 complete lidar points with metrics and derived forest attributes
  ## create dataframe from spatialpointsdataframe
  ## merge with above dataframe of lidar metrics and forest attributes
  ## count number of lidar points in each POLYID
  ## keep only plots with 9 lidar points
  lidar.sample.df <- lidar.sample@data  # starting point is the already "9 plots/polyg"-checked samples but which needs to be rechecked after the merge with the attributes info.
  lidar.metrics.attributes.sample <- merge(lidar.metrics.attributes, lidar.sample.df, by.x = "unique_id", by.y = "Unique_ID")
  dim(lidar.metrics.attributes.sample) ## 45380 points with metrics and attributes out of 50172 total sampled points
  
  npts.sample.plots <- ddply(lidar.metrics.attributes.sample, .(POLY250ID), summarize, npts=length(POLY250ID))
  
  if (paramsGL$unit.size == 75) {
    if ( all(unique(npts.sample.plots$npts) %in% c(0:9)) ) {
      freq.thresh <- 9
    } else {
      stop("Nr. plots per polygon exceeds 9 (with a 75m unit)")  
    }
  } else if (paramsGL$unit.size == 25) {
    if ( all(unique(npts.sample.plots$npts) %in% c(0:1)) ) {
      freq.thresh <- 1
    } else {
      stop("Nr. plots per polygon exceeds 1 (with a 25m unit)")  
    }
  }
  pts9.sample.plots <- subset(npts.sample.plots, npts == freq.thresh)
  nrow(pts9.sample.plots) ## 4372 plots have all 9 lidar points
  nrow(pts9.sample.plots)/nrow(npts.sample.plots) ## 81% of sample lidar plots have all 9 lidar points
  
  ## subset to only sample plots with all 9 lidar points
  pts9.sample.metric.attributes <- merge(lidar.metrics.attributes.sample, pts9.sample.plots, by.x. = "POLY250ID", by.y = "POLY250ID")
  dim(pts9.sample.metric.attributes) ## 39348 points in 4372 plots
  
  
  ############################################################################
  ##remove plots with high within plot heterogeniety using CV of elev_p95 greater than 50% as cutoff
  if (paramsGL$unit.size == 75) {
    cv.elev.p95 <- ddply(pts9.sample.metric.attributes, .(POLY250ID), summarize, cv_elev_p95 = 100 * (sd(elev_p95)/mean(elev_p95)))
  } else if (paramsGL$unit.size == 25) {
    cv.elev.p95 <- ddply(pts9.sample.metric.attributes, .(POLY250ID), summarize, cv_elev_p95 = 100 * elev_cv)
  }
  
  ## histogram of cv_elev_p95
  str <- file.path(wkg_dir, paste(zone, "_hist_cv_elev.pdf", sep=''))
  pdf(str)
    hist(cv.elev.p95$cv_elev_p95, breaks=20)
  dev.off()
  
  ## remove plots with cv_elev_p95 greater than 50% CV of elev_p95 keep only POLY250ID
  sample.POLY250ID <- subset(cv.elev.p95, cv_elev_p95 < params1$CV.thresh, select = POLY250ID) ##4340 combined trainind and validation plots
  str(sample.POLY250ID)
  
  ################################################################################
  ## randomly select 75% 3255 as training plots and 25% 1085 as validation plots
  ## create indicator variable or training versus validation 
  training.POLY250ID <- sample(sample.POLY250ID$POLY250ID, length(as.numeric(sample.POLY250ID$POLY250ID))*params1$trn.pct, replace = FALSE, prob = NULL)
  validation.POLY250ID <- sample.POLY250ID[(!(sample.POLY250ID$POLY250ID %in% training.POLY250ID)),] 
  train.temp <-  data.frame(training.POLY250ID, as.factor(rep("TRAINING", times = length(training.POLY250ID ))))
  valid.temp <-  data.frame(validation.POLY250ID, as.factor(rep("VALIDATION", times = length(validation.POLY250ID ))))
  names(train.temp) <- c("POLY250ID", "TV")
  names(valid.temp) <- c("POLY250ID", "TV")
  train.valid.POLY250ID <- rbind(train.temp , valid.temp)
  
  ############################################################################
  ## create shapefiles for training and validation
  ## polygons for training and validation
  ## polygon centerpoints for training and validation
  ## 9 lidar points within polygons for training and valdiation
  
  ## plot polygons
  poly.training.validation.temp <- merge(lidar.sample.polygons, train.valid.POLY250ID, by.x = "id", by.y = "POLY250ID")
  poly.training.validation.temp.2 <- subset(poly.training.validation.temp, TV == "TRAINING" | TV == "VALIDATION")
  length(poly.training.validation.temp.2)
  names(poly.training.validation.temp.2) <- c("POLY250ID", "TV") 
  writeOGR(poly.training.validation.temp.2 , wkg_dir, paste(zone, "_poly_250m_training_validation", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## center points of plot polygons
  names(lidar.sample.polygons.centerpoints )
  cpt.poly.training.validation.temp <- merge(lidar.sample.polygons.centerpoints, train.valid.POLY250ID, by.x = "POLY250ID", by.y = "POLY250ID")
  cpt.poly.training.validation.temp.2 <- subset(cpt.poly.training.validation.temp, TV == "TRAINING" | TV == "VALIDATION")
  length(cpt.poly.training.validation.temp.2)
  #----------------  MOST IMPORTANT CENTERPOINTS FILE   -----------------------------
  writeOGR(cpt.poly.training.validation.temp.2 , wkg_dir, paste(zone, "_cpt_poly_250m_training_validation", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  #----------------------------------------------------------------------------------
  
  ## 9pts of plot polygons
  ## additional crosswalk required
  ## also some deleting and renaming of duplicate variables 
  pts9.poly.training.validation.temp <- merge(pts9.sample.metric.attributes, train.valid.POLY250ID, by.x = "POLY250ID", by.y = "POLY250ID")
  pts9.poly.training.validation.temp.2 <- subset(pts9.poly.training.validation.temp, TV == "TRAINING" | TV == "VALIDATION")
  pts9.FCID.POLY250ID.xcwalk <- pts9.poly.training.validation.temp.2[,c("POLY250ID","FCID","TV")] 
  pts9.poly.training.validation.temp.3 <- merge(lidar.sample, pts9.FCID.POLY250ID.xcwalk, by.x = "FCID", by.y = "FCID")
  pts9.poly.training.validation.temp.4 <- subset(pts9.poly.training.validation.temp.3, TV == "TRAINING" | TV == "VALIDATION")
#   pts9.poly.training.validation.temp.5 <- pts9.poly.training.validation.temp.4[,-(15)]  # without column nr 15, WHY 15 COLUMNS??????
#   names(pts9.poly.training.validation.temp.5)[14]<- "POLY250ID"
  pts9.poly.training.validation.temp.5 <- subset(pts9.poly.training.validation.temp.4, select=-c(POLY250ID.y))  # without column nr 15, WHY 15 COLUMNS?????? I guess it is to remove the duplicate of POLY250ID.x and .y
  names(pts9.poly.training.validation.temp.5)[6] <- "POLY250ID"
  writeOGR(pts9.poly.training.validation.temp.5 , wkg_dir, paste(zone, "_pts9_poly_250m_training_validation", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  
  ###############################################
  ## get polygon level 9 pt mean values for lidar metrics and derived forest attributes for training and validation
  
  ###############################################################################
  ## create polygon level averages of lidar metrics and derived forest attributes for training and validation
  ## first subset to only polygons used in training and validation
  ## then subset to only POLY250ID and variables to average
  ## then average by POLY250ID
  ## then add indicator variable for training and validation data
  ## finally export as csv
  pts9.metric.attributes.training.validation <-  merge(pts9.sample.metric.attributes, pts9.FCID.POLY250ID.xcwalk, by.x = "FCID", by.y = "FCID")
  pts9.metric.attributes.training.validation.2 <- subset(pts9.metric.attributes.training.validation, TV == "TRAINING" | TV == "VALIDATION")
  pts9.metric.attributes.training.validation.3 <- pts9.metric.attributes.training.validation.2[,c(2,4:87)]
  names(pts9.metric.attributes.training.validation.3)[1] <- "POLY250ID"
  if (paramsGL$unit.size == 75) {
      pts9.mean.metric.attributes.training.validation <- ddply(pts9.metric.attributes.training.validation.3, "POLY250ID", colwise(mean))  # ddply() is a function taking a dataframe (d) as input, doing some operations, and returning another dataframe (d)
  } else if (paramsGL$unit.size == 25) {
      pts9.mean.metric.attributes.training.validation <- pts9.metric.attributes.training.validation.3   # no need to have mean if there is only central plot
  }
  
  merge(pts9.mean.metric.attributes.training.validation, train.valid.POLY250ID, by.x = "POLY250ID", by.y = "POLY250ID")

#   ## clock UTM zone time
#   temp.toc <- proc.time()-temp.tic[3]
#   print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))

}

write.csv(full.df, file = file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", sep = ''))

stopCluster(cl)


# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog1, total elapsed time:",seconds_to_period(toc[3])))