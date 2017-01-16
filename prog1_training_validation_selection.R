## Project Name: NationalImputationForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog1_training_validation_selection.R                          
## Objective: LiDAR training and validation plot selection
##            Merge lidar metrics and derived forest attributes 
##            Remove plots will high variablilty of vegetion height to avoid plots with mixed conditions
##            Generate shapefile of 75 m sqaure polygons of training and validation plots
##            Generate shapefile of center point of 75 m sqaure polygons of training and validation plots
##            Generate shapefile of all 9 points in 75 m sqaure polygons of training and validation plots
##            Calculate mean lidar vegetation metrics and derived forest attributes for all training and validation plots

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - use foreach
# - check no redefinition of paramsGL$zones 
# - check parameters

## SOLVED:
# -V removed "lidar.pts.fc" as samples of the LiDAR transect are already loaded
# -V UTM11S_plot_cloud_metrics_first_returns_gt2m.csv is only 7000 KB in size and Unique_IDs do not match -- with new data from Geordie now it is OK, 7000 Kb bc is a very small transect
# -V where shall I overwrite the new shapefiles from "ForGiona.zip" (use these as there was an issue wrt nr. LiDAR plots >= nr. plots with LiDAR attr. >= nr. plots with forest attr.)? -- In LOP_attributed is OK bc now they have for_sum field
# -V consider doing every operation with poly_filtering and pt_9samples_filtering and only after that take the centerpoint and smaller polygon to produce the final files -- Not done it is much easier to just load "lidar.sample.polygons" as "_poly_sampling". This way "poly_filtering" is not used bc filtering is done in this script with the point shp with 9 plots
# -V check what we keep at this line bc selection based on col indexes: pts9.metric.attributes.training.validation.3 <- pts9.metric.attributes.training.validation.2[,c(2,4:87)] -- adapted to keep Ecozone info and FCID for subsetting to have only centerpoint 
# -V elev_cv is ok for in-plot filtering or we need CV on 3x3 polyg ? -- no it is better to use 3x3 CV, as in Harold's version of the scripts
# -V lidar.metrics.attributes.sample <- merge(lidar.metrics.attributes, lidar.sample.df, by.x = "unique_id", by.y = "Unique_ID"): is empty bc unique_id and Unique_ID do not match -- not the case anymore
# -V TRN/VAL center of polygon shapefiles match with those produced by prog0, the shapefile with all the points written at end  
#   (pts9.poly.training.validation.temp.5 as _pts9_poly_250m_training_validation) does not match!!! -- now it's OK
# -V params1$elev.thresh <- 30 m on elev-p95 is that OK across Canada? let it vary by UTM zone? -- better go higher and set it to 60m not to remove actual trees, also save how many plots are removed in each UTM zone

#### READS/WRITES ------------------------------------------------------------

## READS:
# - "<UTMzone>_plot_cloud_metrics_first_returns_gt2m.csv": csv table with lidar veg metrics database
# - "<UTMzone>_plot_inventory_attributes2.csv": csv table with lidar biomass and inventory attributes database
# - "<UTMzone>_pt_9plots_filtering.shp": (from prog0) point shp with the 9 LiDAR plots to be used for filtering
# - "<UTMzone>_poly_sampling.shp": (from prog0) polygon shp with sampling polygons (3x3 or 1x1) covering plots' squares, after subsetting wrt forest/dist/undist
# - "<UTMzone>_pt_centerpt.shp": (from prog0) point shp with center of hexagons used for sampling, after subsetting wrt forest/dist/undist

## WRITES:
# - "<UTMzone>_poly_training_validation.shp": (for prog2) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)
# - "<UTMzone>_pt_centerpt_training_validation.shp": (for prog2) point shp with centerpoints of selected polygons 
# - "<UTMzone>_pt_9plots_training_validation.shp": point shp with all 3x3 = 9 plots within selected polygons (just to visually check 9 plots configurations, not used anymore later on)
# - "lidar_metrics_mean_training_validation.csv": (for prog3a and prog3) csv table with all observed LiDAR and forest attributes (Y) for selected TRN and VAL samples (3x3 polygons with average plot values or just single plot value) for all UTM zones

#### INIT --------------------------------------------------------------------

print('Prog1, TRN/VAL splitting')

rm(list=ls()) # clear all variables

param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params1 <- list()
params1$elev.thresh <- 60   ## threshold on elev_p95 to eliminate birds and towers
params1$CV.thresh <- 50   ## threshold on coefficient of variation to eliminate plots with too much variability
params1$trn.pct <- 0.75   ## proportion of training samples in the trn/val split

param_file_prog1 = file.path(base_wkg_dir, 'AllUTMzones_params1.Rdata', fsep = .Platform$file.sep) 
save(params1, file = param_file_prog1)

nr.clusters = length(paramsGL$zones)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "spdep",
                      "spatstat",
                      "rgeos",
                      "maptools", 
                      "plyr",
                      "dplyr",   ## to be loaded before foreach to avoid "assertion failed" errors
                      "ggplot2",
                      "data.table",
                      "lubridate", 
                      "doParallel", 
                      "foreach"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time

## check existence and create figures folder for histograms of Elev_CV 
UTMzone.subdir <- file.path(base_figures_dir, "Histograms_UTMzone_level", "Elev_CV", fsep = .Platform$file.sep)
if (! file.exists(UTMzone.subdir)){dir.create(UTMzone.subdir, showWarnings = F, recursive = T)}

## create and write 1st line of log file
log.file <- file.path(base_wkg_dir, "log_prog1.txt", fsep = .Platform$file.sep)
file.create(log.file)
fileConn <- file(log.file, "a")  ## open connection to file and specify we want to append new lines
writeLines(sprintf("Percentage of plots filtered bc elev_p95 > %s m:", params1$elev.thresh), fileConn)  
close(fileConn)

cl <- makeCluster(nr.clusters)
registerDoParallel(cl)
## assigns to full.df the row-wise binding of the last unassigned object (dataframe) of the loop (in this case the one formed with "merge(pts9.mean.me....")
full.df <- foreach (z = 1:length(paramsGL$zones), .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
  ## WHEN NOT USING FOREACH, TESTING PHASE
# for (z in 1:length(paramsGL$zones)) {
  ## WHEN NOT USING FOREACH, TESTING PHASE
  zone <- paramsGL$zones[z]
  wkg_dir <- file.path(base_wkg_dir, zone, fsep = .Platform$file.sep)

#### READ & CHECK FILES -----------------------------------------------------
  
  ## read lidar veg metrics database (metrics from first returns above 2m)
  infile <- file.path(LOP_dir,'LOP_attr','plot_cloud_metrics_first_returns_gt2m',paste(zone,"_plot_cloud_metrics_first_returns_gt2m.csv",sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.1streturns <- read.csv(infile,head=TRUE,sep=",")
  if ( !nrow(lidar.metrics.1streturns) == length(unique(lidar.metrics.1streturns$unique_id)) ) {
    stop(sprintf("Duplicate rows in file %s_plot_cloud_metrics_first_returns_gt2m.csv", zone))
  }

  ## read biomass/inventory attributes database
  infile = file.path(LOP_dir,'LOP_attr','plot_inventory_attributes2',paste(zone,"_plot_inventory_attributes2.csv",sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.forest.attributes <- read.csv(infile,head=TRUE,sep=",")
  if ( !nrow(lidar.metrics.forest.attributes) == length(unique(lidar.metrics.forest.attributes$unique_id)) ) {
    stop(sprintf("Duplicate rows in file %s_plot_inventory_attributes2.csv", zone))
  }
  
  ## read shapefiles
  lidar.sample <- readOGR(dsn=wkg_dir, layer=paste(zone, "_pt_9plots_filtering", sep='')) ## 3x3 lidar plots (contained in the 75m polygons) used for filtering
  lidar.sample.polygons <- readOGR(dsn=wkg_dir, layer=paste(zone,"_poly_sampling", sep='')) ## smallest polyg subset, the one with only 9 plots/polyg
  ## MOST IMPORTANT FILE: CENTERPOINTS FILE
  lidar.sample.polygons.centerpoints <- readOGR(dsn=wkg_dir, layer=paste(zone,"_pt_centerpt", sep='')) ## polygon centerpoints, used instead of lidar.hex.centerpoints to produce shp "UTMXXX_cpt_poly_training_validation"
  ## MOST IMPORTANT FILE: CENTERPOINTS FILE
  
  ## check if coordinate systems are the same between shapefiles
  if ( !all(sapply(list(proj4string(lidar.sample.polygons), proj4string(lidar.sample.polygons.centerpoints)), FUN=identical, proj4string(lidar.sample))) ) {
    stop(sprintf("%s: projections of shapefiles loaded in prog1 do not match.", zone))
  }

#### FILTER WRT elev-p95 (MAX & CV), 9PTS/POLYG AFTER MERGE WITH ATTR ---------------
  
  ## remove lidar plots were elev-p95 is greater than a threshold (indicating birds and towers)
  pct.plots.filtr.elevp95 <- sum(lidar.metrics.1streturns$elev_p95 > params1$elev.thresh)/nrow(lidar.metrics.1streturns)  ## percent of plots filtered due to extreme elev_p95
  fileConn <- file(log.file, "a")  ## open connection to file and specify we want to append new lines
  writeLines(sprintf("%s: %s", zone, pct.plots.filtr.elevp95), fileConn)  ## write to file percentage of plots filtered for each UTM zone to check for weird results by UTM zone
  close(fileConn)
  lidar.metrics.1streturns <- subset(lidar.metrics.1streturns, elev_p95 <= params1$elev.thresh) 
  lidar.metrics.attributes <- merge(lidar.metrics.1streturns, lidar.metrics.forest.attributes, by.x =  "unique_id", by.y = "unique_id")  ## merge lidar metrics and forest attributes by unique_id keeping only the plots with forest attributes
  rm(lidar.metrics.1streturns)    ## remove original two data frames bc they take up too much space
  rm(lidar.metrics.forest.attributes)
  gc()

  ## large number of plots have lidar metrics but no derived forest attributes because only a subset of the plots are in forest condtion class as detemined by sk_mask
  ## so we check how many of the lidar sample plots have 9 complete lidar points with metrics and derived forest attributes
  lidar.sample.df <- lidar.sample@data  ## starting point is the already "9 plots/polyg"-checked samples but which needs to be rechecked after the merge with the attributes info.
  lidar.metrics.attributes.sample <- merge(lidar.metrics.attributes, lidar.sample.df, by.x = "unique_id", by.y = "Unique_ID")     ## merge with above dataframe of lidar metrics and forest attributes
  npts.sample.plots <- ddply(lidar.metrics.attributes.sample, .(POLY250ID), summarize, npts=length(POLY250ID))    ## count number of lidar points in each POLYID
  freq.thresh <- 9    ## keep only plots with 9 lidar points
  pts9.sample.plots <- subset(npts.sample.plots, npts == freq.thresh)     ## subset to only sample plots with all 9 lidar points
  pts9.sample.metric.attributes <- merge(lidar.metrics.attributes.sample, pts9.sample.plots, by.x. = "POLY250ID", by.y = "POLY250ID")
  
  ## remove lidar plots were CV of elev_p95 is greater than 50%
  cv.elev.p95 <- ddply(pts9.sample.metric.attributes, .(POLY250ID), summarize, cv_elev_p95=100*(sd(elev_p95)/mean(elev_p95)))   ## compute CV of elev_p95
  str <- file.path(UTMzone.subdir, paste(zone, "_hist_cv_elev.pdf", sep=''))    ## save the histogram of cv_elev_p95
  pdf(str)
    hist(cv.elev.p95$cv_elev_p95, breaks=20)
  dev.off()
  sample.POLY250ID <- subset(cv.elev.p95, cv_elev_p95 < params1$CV.thresh, select = POLY250ID)   ## keep only plots with cv_elev_p95 lower than 50% and keep only POLY250ID

#### TRN/VAL SPLIT INDEX ------------------------------------------------------------
  
  ## create indicator variable of training/validation (75/25% split)
  set.seed(paramsGL$global.seed)
  training.POLY250ID <- sample(sample.POLY250ID$POLY250ID, length(as.numeric(sample.POLY250ID$POLY250ID))*params1$trn.pct, replace=FALSE, prob=NULL)
  validation.POLY250ID <- sample.POLY250ID[(!(sample.POLY250ID$POLY250ID %in% training.POLY250ID)),] 
  train.temp <- data.frame(training.POLY250ID, as.factor(rep("TRAINING", times=length(training.POLY250ID ))))
  valid.temp <- data.frame(validation.POLY250ID, as.factor(rep("VALIDATION", times=length(validation.POLY250ID ))))
  names(train.temp) <- c("POLY250ID", "TV")
  names(valid.temp) <- c("POLY250ID", "TV")
  train.valid.POLY250ID <- rbind(train.temp , valid.temp)

#### CREATE TRN/VAL SHPS ----------------------------------------------------------

  ## polygons for training and validation
  poly.training.validation.temp <- merge(lidar.sample.polygons, train.valid.POLY250ID, by.x = "id", by.y = "POLY250ID")   ## subsetting with poly.training.validation.temp.2 <- subset(poly.training.validation.temp, TV == "TRAINING")  in Harold's version bc he only built the training set here
  poly.training.validation.temp.2 <- subset(poly.training.validation.temp, TV == "TRAINING" | TV == "VALIDATION")   ## to get rid of NA in TV column (polygons with no that did not pass the filtering by elevation CV)
  names(poly.training.validation.temp.2) <- c("POLY250ID", "TV") 
  writeOGR(poly.training.validation.temp.2, wkg_dir, paste(zone, "_poly_training_validation", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## polygon centerpoints for training and validation
  cpt.poly.training.validation.temp <- merge(lidar.sample.polygons.centerpoints, train.valid.POLY250ID, by.x = "POLY250ID", by.y = "POLY250ID")   ## subsetting with poly.training.validation.temp.2 <- subset(poly.training.validation.temp, TV == "TRAINING")  in Harold's version bc he only built the training set here
  cpt.poly.training.validation.temp.2 <- subset(cpt.poly.training.validation.temp, TV == "TRAINING" | TV == "VALIDATION")
  ## MOST IMPORTANT CENTERPOINTS FILE
  writeOGR(cpt.poly.training.validation.temp.2, wkg_dir, paste(zone, "_pt_centerpt_training_validation", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  ## MOST IMPORTANT CENTERPOINTS FILE
  
  ## 9 lidar points within polygons for training and valdiation
  ## additional crosswalk required, also some deleting and renaming of duplicate variables 
  pts9.poly.training.validation.temp <- merge(pts9.sample.metric.attributes, train.valid.POLY250ID, by.x = "POLY250ID", by.y = "POLY250ID")
  pts9.poly.training.validation.temp.2 <- subset(pts9.poly.training.validation.temp, TV == "TRAINING" | TV == "VALIDATION")
  pts9.FCID.POLY250ID.xcwalk <- pts9.poly.training.validation.temp.2[,c("POLY250ID","FCID","TV")] 
  pts9.poly.training.validation.temp.3 <- merge(lidar.sample, pts9.FCID.POLY250ID.xcwalk, by.x = "FCID", by.y = "FCID")
  pts9.poly.training.validation.temp.4 <- subset(pts9.poly.training.validation.temp.3, TV == "TRAINING" | TV == "VALIDATION")
  pts9.poly.training.validation.temp.5 <- subset(pts9.poly.training.validation.temp.4, select=-c(POLY250ID.y))  ## Harold's version: remove column nr 15 to remove the duplicate of POLY250ID.x and .y
  names(pts9.poly.training.validation.temp.5)[18] <- "POLY250ID"
  writeOGR(pts9.poly.training.validation.temp.5 , wkg_dir, paste(zone, "_pt_9plots_training_validation", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## just to check validity of other layers, never used later on

#### GET TRN/VAL POLYGON VALUES ----------------------------------------------------------

  ## create polygon level values (averages or centerpoint) of lidar metrics and derived forest attributes for training and validation
  pts9.metric.attributes.training.validation <- merge(pts9.sample.metric.attributes, pts9.FCID.POLY250ID.xcwalk, by.x = "FCID", by.y = "FCID")
  pts9.metric.attributes.training.validation.2 <- subset(pts9.metric.attributes.training.validation, TV == "TRAINING" | TV == "VALIDATION")    ## first subset to only polygons used in training and validation
  pts9.metric.attributes.training.validation.3 <- pts9.metric.attributes.training.validation.2[, c(1:2, 100, 4:87)]   ## subset to only POLY250ID and variables to average, also keep FCID (idx=1) and Ecozone column (idx=100)
  names(pts9.metric.attributes.training.validation.3)[2] <- "POLY250ID"
  pts9.metric.attributes.training.validation.3$FCID <- as.character(pts9.metric.attributes.training.validation.3$FCID)  ## avoid unmatching levels if considered as a factor
  pts9.metric.attributes.training.validation.3$POLY250ID <- as.character(pts9.metric.attributes.training.validation.3$POLY250ID)
  if (paramsGL$unit.size == 75) {   ## !!!!!! if we use a unit of 75m we still have to adapt to avoid taking the mean of FCIDs and Ecozone IDs !!!!!!
      pts9.mean.metric.attributes.training.validation <- ddply(pts9.metric.attributes.training.validation.3, "POLY250ID", colwise(mean))  ## average by POLY250ID, ddply() is a function taking a dataframe (d) as input, doing some operations, and returning another dataframe (d)
  } else if (paramsGL$unit.size == 25) {  ## if unit is 25m just take central plot whose FCID is equal to POLY250ID
      pts9.mean.metric.attributes.training.validation <- subset(pts9.metric.attributes.training.validation.3, FCID==POLY250ID)  ## no need to have mean if there is only central plot
  }
  merge(pts9.mean.metric.attributes.training.validation, train.valid.POLY250ID, by.x = "POLY250ID", by.y = "POLY250ID") ## add indicator variable for training and validation data and assign dataframe to full.df

}

## write after sortrows (arrange) based on FCID (for direct comparison with poly_training_validation_exvars_extract.csv)
full.df$FCID <- as.character(full.df$FCID)
write.csv(arrange(full.df, FCID), file = file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", sep = ''))

stopCluster(cl)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog1, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)