########################################################
## Saskatchewan explanatory variable extraction
## extract mean values for all explanatory rasters spatially coincident with sampled lidar plot polygons
## extract is for both training and validation plots
## mean plot level values of explanatory variables is wieghted average based on relative proportion of each pixel within the polygon (aka normalized wieghts)
#######################################################
##----------------------
## TO DO
##----------------------
# - tests to see how extract() works (RAM used, cores, temp directory, etc.)
# - script sequentially copy/pasting files in E from X (only .dat and .hdr)
# - make it general (all zone numbers and S, N, A) then use == UTMzone (my 18) and exist() to check for valid files to copy
# - maybe it is faster to use cellFromXY to get index of pixel (or 4 pixels neighborhood) where plot center is located and query rasters just by indexing (implement weighted average as Piotr did)
# - check extents
# - differences in the rasters read here between Harold's files (eg. "Mosaic_PreChange_persistence_2010.img") and NTEMS ones 
#   (eg."SRef_13S_PreChange_persistence.dat"). Apparently not only in terms of the extent but also in terms of pixel values.
# - check projections to ensure we're sampling proper locations
# - Where do I find Tasseled Cap for all UTM zones?
# - Keep following predictors? 
#     disturbed -- No, info included in "First_Change_Year"
#     Years_Since_Greatest_Change, Years_Since_Last_Change  Years_Since_First_Change -- No, included in "_Change_Year" variables.
#     Aspect, Aspect_cos -- ???
# - Include following NTEMS layers? 
#     SRef_13S_Combined_changes_Multiyear.dat -- No, redundant as this info is used to build all the change layers
#     UTM13S_Change_attribution_complete.dat -- Yes, to add as this is the output of classif of type of changes (dummy vars from it!)
# - xyFromCell() to get coordinates
# -V why weights = FALSE ? -- was wrong in Harolds version: should be set as TRUE to have weighted averages
# -V why not extracting VAL data at the same time?  Employed only at the very end to assess models' accuracy: see "grep -r -n 'TV == "VALIDATION"' *" run in "cd X:/Harolds_orginal_work/SK_nn_forest/sk_docs_code/Rcode"
# - rastertemp directory unused 
# -V script between this and rf_gnn should be the one that takes train and val points and merges them by ecozone -- we don't do that
# - check automatic usage of cores (cluster) in extract()
# - ############# = TODEL
# - check weighted average over 25m unit for all layers



##----------------------
## READS  
##----------------------
# - "<UTMzone>_cpt_poly_250m_training_validation.shp": (from prog1) point shp with centerpoints of selected polygons
# - "<UTMzone>_poly_250m_training_validation.shp": (from prog1) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)


##----------------------
## WRITES
##----------------------
# - "poly_training_validation_exvars_extract.csv": csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots polygons with weighted average pixels values)

#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

rm(list=ls()) # clear all variables

##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.R"
load(param_file)

##----------------------------
## SCRIPT SPECIFIC PARAMETERS
##----------------------------
TC.var.names.long <- list('TC_Brightness', 'TC_Greenness', 'TC_Wetness', 'TC_Angle', 'TC_Distance')
TC.var.names.short <- list('TCB', 'TCG', 'TCW', 'TCA', 'TCD')

cng.var.names.long <- list("PreChange_persistence", "PreChange_magnitude_variation", "PreChange_evolution_rate", 
                        "Change_persistence", "Change_magnitude_variation", "Change_rate",
                        "PostChange_persistence", "PostChange_magnitude_variation", "PostChange_evolution_rate",
                        "Greatest_Change_Year", "Last_Change_Year", "Last_Change_persistence", "First_Change_Year", "First_Change_persistence"
                        )
cng.var.names.short <- list("PreCh_pers", "PreCh_mag", "PreCh_er", 
                        "Ch_pers", "Ch_mag", "Ch_er",
                        "PostCh_pers", "PostCh_mag", "PostCh_er",
                        "GreatCh_yr", "LastCh_yr", "LastCh_pers", "FirstCh_yr", "FirstCh_pers"
                        )

topo.var.names.long <- list("Elevation", "Slope", "Topographic_wetness_index", "Topographic_solar_radiation_index")
topo.var.names.short <- list("Elev", "Slope", "TWI", "TSRI")

params2 <- list()
params2$var.names.long <- c(TC.var.names.long, cng.var.names.long, topo.var.names.long)
params2$var.names.short <- c(TC.var.names.short, cng.var.names.short, topo.var.names.short)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.R', fsep = .Platform$file.sep) 
save(params2, file = param_file_prog2)

nr.clusters = length(paramsGL$zones)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end


##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "landsat",
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
  
  temp.tic <- proc.time() # start clocking time for each UTM zone
  
  zone <- paramsGL$zones[z]
  zone.nr <- substr(zone, 4, nchar(zone))
  print(paste('Prog2, Explanatory variables extraction on' ,zone))   # converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)
  results_dir = file.path(base_results_dir, zone, fsep = .Platform$file.sep)
  
  setwd(wkg_dir)
  
  # ## set basepath
  # basepath<-"F:/SK_nn_forest/sk_data/lidar_plots_database"
  # setwd(basepath)

  ## running large rasters results in the hard drive filling up with temp gribs
  ## deal with it by telling raster package to write temp files in a new defined directory within the working directory
  ## then remove the directory at the end of session
  raster_tmp_dir <- "raster_tmp"
  dir.create(raster_tmp_dir, showWarnings = F, recursive = T)
  rasterOptions(tmpdir = raster_tmp_dir)
  
  ##load polygon centerpoint shapefile to get FCID
  ##load training and validation polygons
  ##subset to only training polygons
  cpt.poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_cpt_poly_250m_training_validation",sep = '')) ## 4340 polygons
  
  poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_poly_250m_training_validation",sep = '')) ## 4340 polygons

    ##extract test area from descending node WRS2
  ##wrs2_descending <-  readOGR(dsn = "F:/SK_nn_forest/sk_data/wrs2", layer = "wrs2_descending") ## 288892 polygons
  ##subset to only wrspr 39022
  ##wrs2_39022 <- subset(wrs2_descending, WRSPR == 39022 )
  ##check proj
  ##proj4string(wrs2_39022)
  ##proj4string(TC)
  ##reproject transform  to match rasters 
  ##wrs2_39022 <- spTransform(wrs2_39022, CRS(proj4string(TC)))
  ##write shapefile
  ##writeOGR(wrs2_39022 , "J:\\sk_data\\wrs2", "wrs2_39022", driver="ESRI Shapefile")
  
  ##--------------------------------
  ##         load rasters
  ##--------------------------------
  ## tassel cap BGW angle and distance (5 names)
  TC <- stack(file.path(TC_dir, paste(zone,"_TC_Mosaic_2010.img",sep = ''),fsep = .Platform$file.sep))
  
  ## disturbance binary mask
  # Disturbed <- raster(file.path(change_dir, "Disturbed_binary_2010.img", fsep = .Platform$file.sep))  # not in Zald RSE 2016
  
  ## change metrics (14 names of files in NTEMS Results/Changes folder)
  ## read raster layers from folder <UTM_zone.nr>\Results\Change_metrics\
  ## change_dir goes at the end as a variable bc if left inside sprintf as a string is not detected by foreach and thus not passed to the cores (error: "object 'change_dir' not found")
  for (var in cng.var.names.long) {
    cmd <- sprintf('%s <- raster(file.path("%s", "UTM_%s", "Results", "Change_metrics", "SRef_%s_%s.dat", fsep = .Platform$file.sep))', var, change_dir, zone.nr, zone.nr, var)  
    eval(parse(text=cmd))
  }
  
  ## topographic metrics (4 names)
  Elevation <- raster(file.path(topo_dir, "DEM_UTM", paste(zone, "_DEM.dat", sep= ''), fsep = .Platform$file.sep))
  Slope <- raster(file.path(topo_dir, "DEM_METRICS", paste(zone, "_TopoMetrics.dat", sep= ''), fsep = .Platform$file.sep), band = 1)  ## band 1 is Slope whereas band 2 and 3 are Aspect and Shaded Relief
  Topographic_wetness_index <- raster(file.path(topo_dir, "TWI", paste(zone, "_TWI.dat", sep= ''), fsep = .Platform$file.sep))
  Topographic_solar_radiation_index <- raster(file.path(topo_dir, "TSRI",  paste(zone, "_TSRI.dat", sep= ''), fsep = .Platform$file.sep))

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
  
  vars <- c(cng.var.names.short)  ## to rename columns of final dataframe
  str <- paste(c(cng.var.names.long), collapse=", ")  ## creates string of variable names to input in stack()
  
  # vars <- c(cng.var.names.short, topo.var.names.short) 
  # str <- paste(c(cng.var.names.long, topo.var.names.long), collapse=", ")
  
  # vars <- params2$var.names.short
  # str <- paste(c("TC", cng.var.names.long, topo.var.names.long), collapse=", ")
  
  cmd <- sprintf("exvars.stack <- stack(%s)", str)
  eval(parse(text=cmd))
  
  ## check number and names of layers in raster brick should be 26
  ## rename layers in raster brick
  class(exvars.stack)
  names(exvars.stack)
#   vars <- c('TCB', 'TCG', 'TCW', 'TCA', 'TCD', 'Disturbed', 'PreChange_persistence', 'PreChange_magnitude', 'PreChange_evolution', 'PostChange_persistence',
#                            'PostChange_magnitude', 'PostChange_evolution', 'Greatest_Change_Year', 'Last_Change_Year', 'Last_Change_persistence', 
#                            'First_Change_Year', 'First_Change_persistence', 'Change_rate', 'Change_persistence', 'Change_magnitude', 
#                            'Years_Since_First_Change', 'Years_Since_Greatest_Change', 'Years_Since_Last_Change',
#                            'Elevation', 'Aspect', 'Aspect_cos', 'Slope', 'TWI', 'TSRI')

  
  
  names(exvars.stack) <- vars
  
  
  ##extract mean of exvars weighted by nomralized proportion contribution of each cell in polygon
  ##weights within each polygon sum to 1
  exvars.extract <- extract(exvars.stack, poly.training.validation, fun = mean, na.rm = FALSE, weights = TRUE, normalizeWeights = TRUE,
                          cellnumbers = FALSE, small = FALSE, df = TRUE, factors = FALSE, sp = FALSE)  # extract is a function that actually accesses the values of the rasters, MEMORY ISSUES?
  
  
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
  
#   ##check the change metrics seperately for disturbed and undisturbed polygons
#   exvars.extract.disturb <- subset(exvars.extract, Disturbed == 1)
#   exvars.extract.undisturb <- subset(exvars.extract, Disturbed == 0)
#   
#   par(mfrow=c(1,2))
#   hist(exvars.extract.disturb$PreChange_persistence)
#   hist(exvars.extract.undisturb$PreChange_persistence)
#   
#   hist(exvars.extract.disturb$PreChange_magnitude)
#   hist(exvars.extract.undisturb$PreChange_magnitude)
#   
#   hist(exvars.extract.disturb$PreChange_evolution)
#   hist(exvars.extract.undisturb$PreChange_evolution)
#   
#   hist(exvars.extract.disturb$PostChange_persistence)
#   hist(exvars.extract.undisturb$PostChange_persistence)
#   
#   hist(exvars.extract.disturb$PostChange_magnitude)
#   hist(exvars.extract.undisturb$PostChange_magnitude)
#   
#   hist(exvars.extract.disturb$PostChange_evolution)
#   hist(exvars.extract.undisturb$PostChange_evolution)
#   
#   hist(exvars.extract.disturb$Greatest_Change_Year)
#   hist(exvars.extract.undisturb$Greatest_Change_Year)
#   
#   hist(exvars.extract.disturb$Last_Change_Year)
#   hist(exvars.extract.undisturb$Last_Change_Year)
#   
#   hist(exvars.extract.disturb$Last_Change_persistence)
#   hist(exvars.extract.undisturb$Last_Change_persistence)
#   
#   hist(exvars.extract.disturb$First_Change_Year)
#   hist(exvars.extract.undisturb$First_Change_Year)
#   
#   hist(exvars.extract.disturb$First_Change_persistence)
#   hist(exvars.extract.undisturb$First_Change_persistence)
#   
#   hist(exvars.extract.disturb$Change_rate)
#   hist(exvars.extract.undisturb$Change_rate)
#   
#   hist(exvars.extract.disturb$Change_persistence)
#   hist(exvars.extract.undisturb$Change_persistence)
#   
#   hist(exvars.extract.disturb$Change_magnitude)
#   hist(exvars.extract.undisturb$Change_magnitude)
#   
#   hist(exvars.extract.disturb$Years_Since_First_Change)
#   hist(exvars.extract.undisturb$Years_Since_First_Change)
#   
#   hist(exvars.extract.disturb$Years_Since_Greatest_Change)
#   hist(exvars.extract.undisturb$Years_Since_Greatest_Change)
#   
#   hist(exvars.extract.disturb$Years_Since_Last_Change)
#   hist(exvars.extract.undisturb$Years_Since_Last_Change)
  
  ##delete temp folder with temp rasters
  unlink(raster_tmp_dir, recursive = T, force = T)
  
  ##bind training polygon information with extracted explanatory variables
  FCID <- cpt.poly.training.validation@data[,c("FCID","POLY250ID","TV")]
  cbind(FCID, exvars.extract)

  
  # clock UTM zone time
#   temp.toc <- proc.time()-temp.tic[3]
#   print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))
  
}

## export full.df (former "poly.training.validation.exvars.extract") as csv
write.csv(full.df, file = file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", sep = ''))

# stopCluster(cl)

# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog2, total elapsed time:",seconds_to_period(toc[3])))
  
  

