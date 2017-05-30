## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog1a_volume_biomass_check_recompute.R                             
## Objective: check the issue with incoherence of Vol and Biomass values between Boreal/Non-boreal ecozones. 
## Wrong equations have been used and this script recomputes the correct values on the transect for these 2 attributes.

#### TO DO -------------------------------------------------------------------

## STILL TO DO:

# Prior to actual run:
# - 

## SOLVED:


#### INIT --------------------------------------------------------------------

print('Prog1a: Volume and biomass check and recompute')

rm(list=ls()) ## clear all variables

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")

#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params1a <- list()

params1a$csv.dir.name.forest.attr.OK <- 'plot_inventory_attributes2'
params1a$csv.dir.name.forest.attr.KO <- 'plot_inventory_attributes2_WRONG_VOL_BIOM'
params1a$outfile.dir <- "D:/Research/ANALYSES/NationalMappingForestAttributes/VolumeCheck"
params1a$outfile.name <- "RMSEs_Bater.rds"

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

csv.dir <- file.path(LOP_dir, 'LOP_attr', fsep = .Platform$file.sep)
csv.file.name.lidar.attr <- "_plot_cloud_metrics_first_returns_gt2m.csv"
csv.file.name.forest.attr <- "_plot_inventory_attributes2.csv"

nr.clusters <- length(paramsGL$zones) 

csv.dir.forest.attr.OK <- file.path(csv.dir, params1a$csv.dir.name.forest.attr.OK, fsep=.Platform$file.sep)
if (! file.exists(csv.dir.forest.attr.OK)){dir.create(csv.dir.forest.attr.OK, showWarnings = F, recursive = T)}

cl <- makeCluster(nr.clusters)
registerDoParallel(cl)
# assigns to full.df the row-wise binding of the last unassigned object (dataframe) of the loop (in this case the one formed with "merge(pts9.mean.me....")
RMSEs <- foreach (z = 1:length(paramsGL$zones), .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
  ## WHEN NOT USING FOREACH, TESTING PHASE
# for (z in 1:length(paramsGL$zones)) {
  ## WHEN NOT USING FOREACH, TESTING PHASE
  
  zone <- paramsGL$zones[z]

#### READ & CHECK FILES -----------------------------------------------------
  
  ## read lidar veg metrics database (metrics from first returns above 2m)
  infile <- file.path(csv.dir, 'plot_cloud_metrics_first_returns_gt2m', paste(zone, csv.file.name.lidar.attr, sep = ''),fsep = .Platform$file.sep)
  lidar.metrics.1streturns <- as.data.table(read.csv(infile,head=TRUE,sep=","))
  setkey(lidar.metrics.1streturns, unique_id)
  
  ## read biomass/inventory attributes database
  infile <- file.path(csv.dir, params1a$csv.dir.name.forest.attr.KO, paste(zone, csv.file.name.forest.attr, sep = ''), fsep=.Platform$file.sep)
  lidar.metrics.forest.attributes <- fread(infile, header=TRUE, sep=",")
  id_order <- copy(lidar.metrics.forest.attributes$unique_id)
  setkey(lidar.metrics.forest.attributes, unique_id)
  
  metrics <- merge(lidar.metrics.1streturns, lidar.metrics.forest.attributes, all=FALSE)
  
  res <- apply_regr_Bater(metrics$elev_mean, metrics$elev_cv, metrics$percentage_first_returns_above_2m)
  
  metrics$gross_stem_volume_OK <- res$gross_stem_volume_OK
  metrics$total_biomass_OK <- res$total_biomass_OK
  setkey(metrics, unique_id)
  
  ## change volume values with new, OK values
  lidar.metrics.forest.attributes <- merge(lidar.metrics.forest.attributes, metrics[, c("unique_id", "gross_stem_volume_OK", "total_biomass_OK"), with=FALSE], all=FALSE)
  lidar.metrics.forest.attributes$gross_stem_volume <- lidar.metrics.forest.attributes$gross_stem_volume_OK
  lidar.metrics.forest.attributes[, gross_stem_volume_OK:=NULL]
  lidar.metrics.forest.attributes$total_biomass <- lidar.metrics.forest.attributes$total_biomass_OK
  lidar.metrics.forest.attributes[, total_biomass_OK:=NULL]
  
  ## reorder rows to be same as when read
  lidar.metrics.forest.attributes <- lidar.metrics.forest.attributes[match(id_order, lidar.metrics.forest.attributes$unique_id), ]
  
  ## save table to CSV into a new, OK folder params1a$csv.dir.name.forest.attr.OK
  outfile <- file.path(csv.dir.forest.attr.OK, paste(zone, csv.file.name.forest.attr, sep = ''), fsep=.Platform$file.sep)
  fwrite(lidar.metrics.forest.attributes, outfile)

  c(sqrt(sum((res$gross_stem_volume - metrics$gross_stem_volume)**2)/nrow(metrics)), 
    sqrt(sum((res$gross_stem_volume_OK - metrics$gross_stem_volume)**2)/nrow(metrics)), 
    sqrt(sum((res$total_biomass - metrics$total_biomass)**2)/nrow(metrics)),
    sqrt(sum((res$total_biomass_OK - metrics$total_biomass)**2)/nrow(metrics)))

}

stopCluster(cl)

rownames(RMSEs) <- paramsGL$zones
colnames(RMSEs) <- c("gross_stem_volume", "gross_stem_volume_OK", "total_biomass", "total_biomass_OK")

outfile <- file.path(params1a$outfile.dir, params1a$outfile.name, fsep=.Platform$file.sep)
saveRDS(RMSEs, outfile)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog1a, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)


