########################################################
## Explanatory variable extraction
## extract mean values for all explanatory rasters spatially coincident with sampled lidar plot polygons for both training and validation plots
## mean plot level values of explanatory variables is weighted average based on relative proportion of each pixel within the polygon (aka normalized wieghts)
####################################################### 
##----------------------
## TO DO
##----------------------

## STILL TO DO
# - convert "Year of " variables to "Years since" to avoid weird averages
# - adapt behavior of params2$pred.names.lg and params2$pred.names.sh wrt list vs vector
# -TXOMIN: check extents (Topo by Geordie vs rest) -- Geordie has a script to change and resave all NTEMS layers at the new narrow extent (less overlap) with the same folder structure
# - check projections to ensure we're sampling proper locations
# - check var names wrt final df (35 vs 32 cols)
# - how to ignore unclassified in Change_attribution? 
# - hours of daylight from lat long -- only a function of Lat, annual averages out, length of day at Aug 1st is better
# -DOUG: climate variables at 1km (finer res?)
# - differences in the rasters read here between Harold's files (eg. "Mosaic_PreChange_persistence_2010.img") and NTEMS ones 
#   (eg."SRef_13S_PreChange_persistence.dat"). Apparently not only in terms of the extent but also in terms of pixel values.
# - Keep following predictors? 
#     disturbed -- No, info included in "First_Change_Year"
#     Years_Since_Greatest_Change, Years_Since_Last_Change  Years_Since_First_Change -- No, included in "_Change_Year" variables.
#     Aspect, Aspect_cos -- ???
# - ############# = TODEL
# - check weighted average over 25m unit for all layers
# - check metric for disturbed and undisturbed polygons... to do (all these histograms that are now removed)?
# - TSRI different values the two different layers
# - rastertemp directory is used?
# - !!! uncomment unlink() !!!
# - !!! uncomment params2$var.names.short !!!
# - prior to actual run, delete all UTMzones folders in Landsat_dir E:\NTEMS otherwise copy_paste_UTMdata() will not be run


## SOLVED
# -V why weights = FALSE ? -- was wrong in Harolds version: should be set as TRUE to have weighted averages
# -V why not extracting VAL data at the same time? Employed only at the very end to assess models' accuracy: see "grep -r -n 'TV == "VALIDATION"' *" run in "cd X:/Harolds_orginal_work/SK_nn_forest/sk_docs_code/Rcode"
# -V script between this and rf_gnn should be the one that takes train and val points and merges them by ecozone -- we don't do that
# -V tests to see how extract() works (RAM used, cores, temp directory, etc.) -- can be implemented in parallel with foreach of blocks of poly.training.validation (19 mins vs 6h30 !!! still to see with all layers but...)
# -V check integration of results of the extraction in full.table, test on zone 13S and 20 since now it works without duplicates -- works fine
# -V script sequentially copy/pasting files in E from X (only .dat and .hdr) -- function copy_paste_UTMdata()
# -V make it general (all zone numbers and S, N, A) then use == UTMzone (my 18) and exist() to check for valid files to copy -- not needed to built TRN set, just in test
# -V maybe it is faster to use cellFromXY to get index of pixel (or 4 pixels neighborhood) where plot center is located and query rasters just by indexing (implement weighted average as Piotr did) -- not needed now
# -V Where do I find Tasseled Cap for all UTM zones? -- To compute it myself from coefficients
# -V xyFromCell() or coordinates() to get coordinates, spTransform() to get lat/long -- coordinates() is the right function
# -V TC transform (coeff for TM vs ETM+ sensors available in 2010?, DN in 16 bits to reflect in 0-1?) -- coefficients corresponding to ETM+, use 'output <- all.bands %*% tc.coef', no need to transf into 0-1 as we don't look at physical meaning, atan2() for TC Angleto automatically select quadrants
# -V Include following NTEMS layers? 
#     SRef_13S_Combined_changes_Multiyear.dat -- No, redundant as this info is used to build all the change layers
#     UTM13S_Change_attribution_complete.dat -- Yes, to add as this is the output of classif of type of changes (no need to create dummy vars from it as both randomForests() and yaimpute+RF can handle categorical inputs)
# -V check automatic usage of cores (cluster) in extract() -- uses 1 core only if not called within a foreach (multiple parallel extract can be run efficiently this way)
# -V how to aggregate change classes (unclass)?  -- leave them like that with 10 classes (and separated 'unclassified' classes), we can always merge them later on
# -V use 500m hexag? -- let's keep 250m and try to model things with a RF on 100'000+ points (parallel!)
# -V no point in saving a dataset (table) for mapping, do a script creating variables on the fly for prediction/mapping -- yes, thats the weay to do it, consider saving a TC components file (5 layers) with the new narrow extent
# -V add ArcGIS/ArcPy AreaSolarRadiation() ? -- no, super-slow to compute


##----------------------
## READS
##----------------------
# - "<UTMzone>_cpt_poly_250m_training_validation.shp": (from prog1) point shp with centerpoints of selected polygons
# - "<UTMzone>_poly_250m_training_validation.shp": (from prog1) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)


##----------------------
## WRITES
##----------------------
# - "poly_training_validation_exvars_extract.csv": (for prog3a and prog3) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

rm(list=ls()) # clear all variables

##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

source("Functions_NatImp.R")

##----------------------------
## SCRIPT SPECIFIC PARAMETERS
##----------------------------

## Prior to actual run, delete all UTMzones folders in Landsat_dir E:\NTEMS otherwise copy_paste_UTMdata() will not be run
params2 <- list()

bands.names.lg <- list('band1', 'band2', 'band3', 'band4', 'band5', 'band7')
bands.names.sh <- list('b1', 'b2', 'b3', 'b4', 'b5', 'b7')

TC.names.lg <- list('TC_Brightness', 'TC_Greenness', 'TC_Wetness', 'TC_Angle', 'TC_Distance')  
TC.names.sh <- list('TCB', 'TCG', 'TCW', 'TCA', 'TCD')

cng.names.lg <- list("PreChange_persistence", "PreChange_magnitude_variation", "PreChange_evolution_rate", 
                        "Change_persistence", "Change_magnitude_variation", "Change_rate",
                        "PostChange_persistence", "PostChange_magnitude_variation", "PostChange_evolution_rate",
                        "Greatest_Change_Year", "Last_Change_Year", "Last_Change_persistence", "First_Change_Year", "First_Change_persistence"
                        )
cng.names.sh <- list("PreCh_pers", "PreCh_mag", "PreCh_er", 
                        "Ch_pers", "Ch_mag", "Ch_er",
                        "PostCh_pers", "PostCh_mag", "PostCh_er",
                        "GreatCh_yr", "LastCh_yr", "LastCh_pers", "FirstCh_yr", "FirstCh_pers"
                        )

cngattr.names.lg <- list("Change_attribution")
cngattr.names.sh <- list("Ch_attr")

topo.names.lg <- list("Elevation", "Slope", "Topographic_wetness_index", "Topographic_solar_radiation_index")
topo.names.sh <- list("Elev", "Slope", "TWI", "TSRI")

trends.names.lg <- list("Longitude", "Latitude", "Hours_of_daylight")
trends.names.sh <- list("Long", "Lat", "Daylight")

params2$pred.names.lg <- list(bands=bands.names.lg, TC=TC.names.lg, cng=cng.names.lg, cngattr=cngattr.names.lg, trends=trends.names.lg)
params2$pred.names.sh <- list(bands=bands.names.sh, TC=TC.names.sh, cng=cng.names.sh, cngattr=cngattr.names.sh, trends=trends.names.sh)


params2$jul.day <- 213  # compute length of daylight for August 1st

## USE THIS TILL ALL LAYERS ARE READY WITH SAME EXTENT
params2$var.names.long <- c(bands.names.lg, TC.names.lg, cng.names.lg, cngattr.names.lg, trends.names.lg)
params2$var.names.short <- c(bands.names.sh, TC.names.sh, cng.names.sh, cngattr.names.sh, trends.names.sh)

##-------------- TO UNCOMMENT ----------------------
# params2$var.names.long <- c(bands.names.lg, TC.names.lg, cng.names.lg, cngattr.names.lg, topo.names.lg, trends.names.lg)
# params2$var.names.short <- c(bands.names.sh, TC.names.sh, cng.names.sh, cngattr.names.sh, topo.names.sh, trends.names.sh)
##-------------- TO UNCOMMENT ----------------------

params2$disk.names <- list("//2234b/j", "//Frst-cdw-2231j/j", "//Frst-cdw-2231j/h", "//Frst-cdw-2231j/f")

## not used so far as nr available cores defines block size
# params2$sizeBlocks <- 200   

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
save(params2, file = param_file_prog2)


##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "landsat",
                      "lattice",
                      "lubridate", 
                      "doParallel", 
                      "foreach", 
                      "data.table"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

tic <- proc.time() # start clocking global time

copy.time.UTMzone <- vector("list", length(paramsGL$zones)) 
full.table <- vector("list", length(paramsGL$zones)) 
for (z in 1:length(paramsGL$zones)) {
  
  temp.tic <- proc.time() # start clocking time for each UTM zone
  
  zone <- paramsGL$zones[z]
  zone.nr <- substr(zone, 4, nchar(zone))
  print(paste('Prog2, Explanatory variables extraction on' ,zone))   # converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)
  results_dir = file.path(base_results_dir, zone, fsep = .Platform$file.sep)
  
  # setwd(wkg_dir)
  
  ## call copy_paste_UTMdata() only if no folder named like that exists in Landsat_dir
  cmd <- sprintf('folder.exists.bool <- file.exists(file.path("%s", "UTM_%s", fsep = .Platform$file.sep))', Landsat_dir, zone.nr)
  eval(parse(text=cmd))
  if ( !folder.exists.bool ) {
    log.copy.paste <- copy_paste_UTMdata(params2$disk.names, zone, cng.names.lg, Landsat_dir)
    if (length(log.copy.paste$list.missing.files)>0) {
      stop(sprintf("Missing files for zone %s", zone))
    }
    copy.time.UTMzone[z] <- paste(zone, "copy time:", seconds_to_period(log.copy.paste$toc[3]))
  }
  
  ## running large rasters results in the hard drive filling up with temp gribs
  ## deal with it by telling raster package to write temp files in a new defined directory within the working directory
  ## then remove the directory at the end of session
  raster_tmp_dir <- "raster_tmp"
  dir.create(raster_tmp_dir, showWarnings = F, recursive = T)  ## creates temp directory in BAP_Imputation_working/wkg/<UTMzone>
  rasterOptions(tmpdir = raster_tmp_dir)
  
  ## load polygon centerpoint shapefile to get FCID
  ## load training and validation polygons
  cpt.poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_cpt_poly_250m_training_validation",sep = '')) ## 4340 polygons
  
  poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_poly_250m_training_validation",sep = '')) ## 4340 polygons

########### TO DEL ###########
#   extract test area from descending node WRS2
#   wrs2_descending <-  readOGR(dsn = "F:/SK_nn_forest/sk_data/wrs2", layer = "wrs2_descending") ## 288892 polygons
#   subset to only wrspr 39022
#   wrs2_39022 <- subset(wrs2_descending, WRSPR == 39022 )
#   check proj
#   proj4string(wrs2_39022)
#   proj4string(TC)
#   reproject transform  to match rasters 
#   wrs2_39022 <- spTransform(wrs2_39022, CRS(proj4string(TC)))
#   write shapefile
#   writeOGR(wrs2_39022 , "J:\\sk_data\\wrs2", "wrs2_39022", driver="ESRI Shapefile")
########### TO DEL ###########
  

  ##--------------------------------
  ##         load rasters
  ##--------------------------------
  
  ########### TO DEL ###########
  ## disturbance binary mask
  # Disturbed <- raster(file.path(change_dir, "Disturbed_binary_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
  ## tassel cap BGW angle and distance (5 names)
  # TC <- stack("X:/ForGiani_WorkingImputation/data/tc_composites/TC_Mosaic_2010.img")
  ########### TO DEL ###########
  
  ## Landsat bands for year 2010 (6 layers), TC layer will be computed later (first for TRN/VAL and later for all Canada when mapping)
  cmd <- sprintf('Bands <- brick(file.path("%s", "UTM_%s", "Results", "proxy_values", "SRef_UTM%s_2010_proxy.dat", fsep = .Platform$file.sep))', Landsat_dir, zone.nr, zone.nr)  
  eval(parse(text=cmd))
  
  ## change metrics (14 layers in NTEMS Results/Changes folder)
  ## read raster layers from folder <UTM_zone.nr>\Results\Change_metrics\
  ## Landsat_dir goes at the end as a variable bc if left inside sprintf as a string is not detected by foreach and thus not passed to the cores (error: "object 'Landsat_dir' not found")
  for (var in cng.names.lg) {
    cmd <- sprintf('%s <- raster(file.path("%s", "UTM_%s", "Results", "Change_metrics", "SRef_%s_%s.dat", fsep = .Platform$file.sep))', var, Landsat_dir, zone.nr, zone.nr, var)  
    eval(parse(text=cmd))
  }
  
  ## change attribution (1 layer)
  cmd <- sprintf('Change_attribution <- raster(file.path("%s", "UTM_%s", "Results", "change_attribution", "%s_Change_attribution_complete.dat", fsep = .Platform$file.sep))', Landsat_dir, zone.nr, zone)  
  eval(parse(text=cmd))
  
  ## topographic metrics (4 layers)
  Elevation <- raster(file.path(topo_dir, "DEM_UTM", paste(zone, "_DEM.dat", sep= ''), fsep = .Platform$file.sep))
  Slope <- raster(file.path(topo_dir, "DEM_METRICS", paste(zone, "_TopoMetrics.dat", sep= ''), fsep = .Platform$file.sep), band = 1)  ## band 1 is Slope whereas band 2 and 3 are Aspect and Shaded Relief
  Topographic_wetness_index <- raster(file.path(topo_dir, "TWI", paste(zone, "_TWI.dat", sep= ''), fsep = .Platform$file.sep))
  Topographic_solar_radiation_index <- raster(file.path(topo_dir, "TSRI",  paste(zone, "_TSRI.dat", sep= ''), fsep = .Platform$file.sep))
  
########### TO DEL ###########
## corresponds to loading the following  
#   PreChange_persistence <- raster(file.path(change_dir, "Mosaic_PreChange_persistence_2010.img", fsep = .Platform$file.sep))
#   PreChange_magnitude <- raster(file.path(change_dir, "Mosaic_PreChange_magnitude_variation_2010.img", fsep = .Platform$file.sep))
#   PreChange_evolution <- raster(file.path(change_dir, "Mosaic_PreChange_evolution_rate_2010.img", fsep = .Platform$file.sep))
#   PostChange_persistence <- raster(file.path(change_dir, "Mosaic_PostChange_persistence_2010.img", fsep = .Platform$file.sep))
#   PostChange_magnitude <- raster(file.path(change_dir, "Mosaic_PostChange_magnitude_variation_2010.img", fsep = .Platform$file.sep))
#   PostChange_evolution <- raster(file.path(change_dir, "Mosaic_PostChange_evolution_rate_2010.img", fsep = .Platform$file.sep))
#   Greatest_Change_Year <- raster(file.path(change_dir, "Mosaic_Greatest_Change_Year_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   Last_Change_Year <- raster(file.path(change_dir, "Mosaic_Last_Change_Year_2010.img", fsep = .Platform$file.sep))   # not in Zald RSE 2016
#   Last_Change_persistence <- raster(file.path(change_dir, "Mosaic_Last_Change_persistence_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   First_Change_Year <- raster(file.path(change_dir, "Mosaic_First_Change_Year_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   First_Change_persistence <- raster(file.path(change_dir, "Mosaic_First_Change_persistence_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   Change_rate <- raster(file.path(change_dir, "Mosaic_Change_rate_2010.img", fsep = .Platform$file.sep))
#   Change_persistence <- raster(file.path(change_dir, "Mosaic_Change_persistence.img", fsep = .Platform$file.sep))
#   Change_magnitude <- raster(file.path(change_dir, "Mosaic_Change_magnitude_variation_2010.img", fsep = .Platform$file.sep))

#   Elevation <- raster(file.path(topo_dir, "CDED.img", fsep = .Platform$file.sep))
#   Aspect <- raster(file.path(topo_dir, "ASPECT_DEG.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   Aspect_cos <- raster(file.path(topo_dir, "ASPECT_COS.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   Slope <- raster(file.path(topo_dir, "SLOPE_DEG.img", fsep = .Platform$file.sep))
#   TWI <- raster(file.path(topo_dir, "TWI.img", fsep = .Platform$file.sep))
#   TSRI <- raster(file.path(topo_dir, "TSRI.img", fsep = .Platform$file.sep))
  
## these files are not in NTEMS folder but convey the same info as "Greatest_Change_Year", "Last_Change_Year" and "First_Change_Year". 
#   Years_Since_Greatest_Change <- raster(file.path(change_dir, "Mosaic_Years_Since_Greatest_Change_2010.img", fsep = .Platform$file.sep))
#   Years_Since_Last_Change <- raster(file.path(change_dir, "Mosaic_Years_Since_Last_Change_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
#   Years_Since_First_Change <- raster(file.path(change_dir, "Mosaic_Years_Since_First_Change_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
  
#   ##make a rasterstack of all potential explantory variables
#   exvars.stack <- stack(TC, Disturbed, PreChange_persistence, PreChange_magnitude, PreChange_evolution, PostChange_persistence,
#                   PostChange_magnitude, PostChange_evolution, Greatest_Change_Year, Last_Change_Year, Last_Change_persistence, 
#                   First_Change_Year, First_Change_persistence, Change_rate, Change_persistence, Change_magnitude, 
#                   Years_Since_First_Change, Years_Since_Greatest_Change, Years_Since_Last_Change,
#                   Elevation, Aspect, Aspect_cos, Slope, TWI, TSRI)
########### TO DEL ###########  
  
  ########### TO DEL ###########  
  #  ## test with padding (pasting the narrower raster into a blank wider one)
  #   unfilled.raster <- raster(extent(Bands), ncol=Bands@ncols, nrow=Bands@nrows, crs=proj4string(Bands))
  #   low.left.coords.unfill <- unfilled.raster@extent[c(1,3)]
  #   low.left.coords.elev <- Elevation@extent[c(1,3)]
  #   
  #   low.left.idx.elev <- round((low.left.coords.elev-low.left.coords.unfill)/30,  digits = 0)
  #   
  #   aaa <- unfilled.raster@data[low.left.idx.elev[2]:low.left.idx.elev[2]+Elevation@nrows, low.left.idx.elev[1]:low.left.idx.elev[1]+Elevation@ncols]
  #   
  #   NOT WORKING CAUSE CANT ACCESS/WRITE DATA IN A RASTER OBJECT
  #   unfilled.raster low.left.corner
  #  
  #   ## test with mask()
  #   cmd <- sprintf('non.overl.mask <- raster(file.path("%s", "UTM_%s", "Non_overlapping_mask", "%s_Non_Overlapping_Mask.dat", fsep = .Platform$file.sep))', Landsat_dir, zone.nr, zone)  
  #   eval(parse(text=cmd))
  #   Bands.masked <- mask(Bands, non.overl.mask, maskvalue=0)
  ########### TO DEL ###########
  
  special.case.names.sh <- c("GreatCh_yr", "LastCh_yr", "FirstCh_yr")
  regular.case.idx <- (!cng.names.sh  %in%  special.case.names.sh)
  
  XXXXXXXXXXXXXXXXXXXXXXXXXXXX
  XXXXXXXXXXXXXXXXXXXXXXXXXXXX
  2 (?) separate extract() to integrate in foreach (rbindlist?)
  "Ch_attr"
  
  
  vars <- c(bands.names.sh, cng.names.sh[regular.case.idx])
  str <- paste(c('Bands', cng.names.lg[regular.case.idx]), collapse=", ")
  
  ## to use till we don't have topo and other NTEMS data matching in extent
  vars <- c(bands.names.sh, cng.names.sh, cngattr.names.sh)
  str <- paste(c('Bands', cng.names.lg, cngattr.names.lg), collapse=", ")
  
#   ## final list to to use when all layers are available
#   vars <- c(bands.names.sh, cng.names.sh, cngattr.names.sh, topo.names.sh)  ## to rename columns of final dataframe
#   str <- paste(c('Bands', cng.names.lg, cngattr.names.lg, topo.names.lg), collapse=", ")  ## creates string of variable names to input in stack()
  
  cmd <- sprintf("exvars.stack <- stack(%s)", str)  ## stack() already checks if projections of all the layers match
  eval(parse(text=cmd))
  names(exvars.stack) <- vars
  
  ## check if coordinate systems are the same between raster stack and shapefiles
  if ( !all(sapply(list(proj4string(poly.training.validation), proj4string(cpt.poly.training.validation)), FUN = identical, proj4string(exvars.stack))) ) {
    stop(sprintf("%s: projections of raster stack and shapefiles loaded in prog2 do not match.", zone))
  }
  
########### TO DEL ###########  
#   vars <- c('TCB', 'TCG', 'TCW', 'TCA', 'TCD', 'Disturbed', 'PreChange_persistence', 'PreChange_magnitude', 'PreChange_evolution', 'PostChange_persistence',
#                            'PostChange_magnitude', 'PostChange_evolution', 'Greatest_Change_Year', 'Last_Change_Year', 'Last_Change_persistence', 
#                            'First_Change_Year', 'First_Change_persistence', 'Change_rate', 'Change_persistence', 'Change_magnitude', 
#                            'Years_Since_First_Change', 'Years_Since_Greatest_Change', 'Years_Since_Last_Change',
#                            'Elevation', 'Aspect', 'Aspect_cos', 'Slope', 'TWI', 'TSRI')
########### TO DEL ###########

  nblocks = min(detectCores(), nrow(poly.training.validation@data))
  sizeBlocks <- ceiling(nrow(poly.training.validation@data)/nblocks)
  nr.clusters <- nblocks
  cl <- makeCluster(nr.clusters)
  registerDoParallel(cl)
  getDoParWorkers()
  exvars.extract <- foreach (blck = 1:nblocks, .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
  # for (blck in 1:nblocks) {
    indPart = seq(from=(blck-1)*sizeBlocks+1, to=min(blck*sizeBlocks, nrow(poly.training.validation@data)), by=1)
    ## extract mean of exvars weighted by nomralized proportion contribution of each cell in polygon
    ## weights within each polygon sum to 1
    ## automatically assigns an ID starting at 1 block of samples assigned to the cores (to be removed later on bc meaningless)
    extract(exvars.stack, poly.training.validation[indPart, ], fun = mean, na.rm = FALSE, weights = TRUE, normalizeWeights = TRUE,
                          cellnumbers = FALSE, small = FALSE, df = TRUE, factors = FALSE, sp = FALSE)  # extract is a function that actually accesses the values of the rasters, MEMORY ISSUES?
  }
  
  stopCluster(cl)
  
#------------ TO UNCOMMENT -----------------  
#   cmd <- sprintf('unlink(file.path("%s", "UTM_%s", fsep = .Platform$file.sep), recursive = T, force = T)', Landsat_dir, zone.nr)
#   eval(parse(text=cmd))
#------------ TO UNCOMMENT -----------------  
  

########### TO DEL ###########  
#   for (var in vars){
#     cmd=sprintf('hist(exvars.extract$%s)', var)
#     eval(parse(text=cmd))
#   }
#   ##check data with histograms
#   hist(exvars.extract$TCB)
#   hist(exvars.extract$TCG)
#   hist(exvars.extract$TCW)
#   hist(exvars.extract$TCA)
#   hist(exvars.extract$TCD)
#   hist(exvars.extract$Elevation)
#   hist(exvars.extract$Aspect)
#   hist(exvars.extract$Aspect_cos)
#   hist(exvars.extract$Slope)
#   hist(exvars.extract$TWI)
#   hist(exvars.extract$TSRI)
########### TO DEL ###########  
  
  
  ##delete temp folder with temp rasters
  unlink(raster_tmp_dir, recursive = T, force = T)
  
  ## coordinates of plot polygons
  long.lat <- as.data.frame(coordinates(spTransform(cpt.poly.training.validation, CRS("+proj=longlat"))))
  daylight <- lat_2_daylight(long.lat[,2], params2$jul.day)
  trends <- as.data.frame(cbind(long.lat, daylight))
  colnames(trends) <- trends.names.sh
  
  UTMzone.ID <- as.data.frame(rep(zone.nr, nrow(cpt.poly.training.validation)))
  colnames(UTMzone.ID) <- "UTMzone"
  
  ##bind training polygon information with UTMzone ID, extracted explanatory variables and trends (coordinates, etc.)
  FCID <- cpt.poly.training.validation@data[,c("FCID","POLY250ID","TV")]
  full.table[[z]] <- cbind(FCID, UTMzone.ID, exvars.extract, trends)
  
  ## clock UTM zone time
  temp.toc <- proc.time()-temp.tic[3]
  print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))
  
}

full.df <- as.data.frame(rbindlist(full.table))

all.bands <- as.matrix(full.df[, c("b1", "b2", "b3", "b4", "b5", "b7")])
TC.list <- bands_2_TC(all.bands)
TC.df <- cbind(TC.list$TCs, TC.list$TCang, TC.list$TCdist)
colnames(TC.df) <- TC.names.sh

## add the 5 new TC components
full.df <- cbind(full.df, TC.df)

## reorder columns to match predefined order and remove column ID
full.df <- full.df[, c("FCID", "POLY250ID", "TV", "UTMzone", unlist(params2$var.names.short))] 

## export full.df (former "poly.training.validation.exvars.extract") as csv
write.csv(full.df, file = file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", sep = ''))

## clock global time
toc <- proc.time()-tic[3]
end.message <- paste("Prog2, total elapsed time:",seconds_to_period(toc[3]))
print(end.message)

## write log file
unlink(file.path(base_wkg_dir, "log_prog2.txt", sep = ''), recursive = T, force = T)
list.to.print <- append(copy.time.UTMzone, end.message)
lapply(list.to.print, write, file.path(base_wkg_dir, "log_prog2.txt", sep = ''), append=TRUE)





