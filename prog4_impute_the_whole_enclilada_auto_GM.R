########################################################
## Saskatchewan RF Imputation Mapping 2010
## Mapping only a subset of SK
## using randomforest model with 500 trees
## single nearest neighbor (k=1) imputation mapping
## 2010 map year
## response variables from lidar plots
## explanatory variables from landsat spectral indices landsat change metrics and topographic indices
## explanatory variables extracted for lidar plots
#######################################################
##----------------------
## TO DO
##----------------------
# - in masking step change name of mask to use again the complete masks and not TRIAL ones
# - check correspondence between variable names in raster (prediction grid) and Xrefs (training set) in model file
# - check correspondence between size of unit (75x75 m for training VS 30x30 m pixel in test)
# - fix units of variables: [kg/ha], etc.
# - getValues() or getValuesBlock() to access reshaped raster values (rows = pixels, columns = bands)


##----------------------
## MAIN OUTPUT
##----------------------
# - script merging "impute_the_whole_encilada.R" and "attaching_plot_data.R", so it writes both the ID raster (result of imputation) and a raster with the forest attribute of interest

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




nr.clusters = length(paramsGL$zones)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end


##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "randomForest",
                      "rgl",
                      "rasterVis",
                      "yaImpute",
                      "snow",
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

# cl <- makeCluster(nr.clusters)
# registerDoParallel(cl)
# getDoParWorkers()
# foreach (z = 1:length(paramsGL$zones), .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
for (z in 1:length(paramsGL$zones)) {
  
  temp.tic <- proc.time() # start clocking time for each UTM zone
  
  zone = paramsGL$zones[z]
  print(paste('Prog4, prediction on' ,zone))   # converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)
  results_dir = file.path(base_results_dir, zone, fsep = .Platform$file.sep)
  
  setwd(wkg_dir)
  
#   ## set basepath
#   basepath<-"G:/ForGiani_WorkingImputation/data"
#   setwd(basepath)
  
  #pnt_dir = file.path(LOP_dir,'LOP_attributed', fsep = .Platform$file.sep) 
  
  
  ## running large rasters results in the hard drive filling up with temp gribs
  ## deal with it by telling raster package to write temp files in a new defined directory within the working directory
  ## then remove the directory at the end of session
  raster_tmp_dir <- "raster_tmp"
  dir.create(raster_tmp_dir, showWarnings = F, recursive = T)
  rasterOptions(tmpdir = raster_tmp_dir)
  
  ##load Rdata rf model and plot data attachment objects from RData file
  # load("G:/ForGiani_WorkingImputation/data/yai_rf_models/yai_rf_gnn_sk.RData") # directory on X:/ disk copied to the following for local use
  load("D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/REFERENCE MODELS/yai_rf_models/yai_rf_gnn_sk_newNames.RData")  # reference model of Zald in local directory
  # load(file.path(wkg_dir, "yai_rf_gnn_sk.RData", fsep = .Platform$file.sep))  # model built in previous step (prog3_rf_gnn...) 
  
    ##load boundary cropping shapefile
#   bnd <- readOGR(dsn = "G:/ForGiani_WorkingImputation/data/masks", layer = "skbnd")
#   bnd <- raster("G:/ForGiani_WorkingImputation/data/masks/skbnd.img")
#   bnd <- readOGR(dsn = file.path(data_dir, "masks", fsep = .Platform$file.sep), layer = "skbnd")
#   bnd <- raster(file.path(data_dir, "masks", "skbnd.img", fsep = .Platform$file.sep))
  # bnd <- readOGR(dsn = file.path(data_dir, "masks", fsep = .Platform$file.sep), layer = "skbnd_TRIAL")  # small subregion to try the script
  bnd <- raster(file.path(data_dir, "masks", "skbnd_TRIAL.img", fsep = .Platform$file.sep))
  
  ##load binary forest mask
#   sk.fc.mask <- raster("G:/ForGiani_WorkingImputation/data/masks/sk_forest_mask.img")
#   fc.mask <- raster(file.path(data_dir, "masks", "sk_forest_mask.img", fsep = .Platform$file.sep))
  fc.mask <- raster(file.path(data_dir, "masks", "sk_forest_mask_TRIAL.img", fsep = .Platform$file.sep))
  
## TO DEL  
#   ##load TCs TCB TCG TCW TCA TCD
#   TC <- stack("G:/ForGiani_WorkingImputation/data/tc_composites/TC_Mosaic_2010.img")
#   ##load change metrics
#   PreChange_persistence <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_PreChange_persistence_2010.img")
#   PreChange_magnitude <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_PreChange_magnitude_variation_2010.img")
#   PreChange_evolution <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_PreChange_evolution_rate_2010.img")
#   PostChange_persistence <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_PostChange_persistence_2010.img")
#   PostChange_magnitude <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_PostChange_magnitude_variation_2010.img")
#   PostChange_evolution <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_PostChange_evolution_rate_2010.img")
#   Years_Since_Greatest_Change <- ("G:/ForGiani_WorkingImputation/data/change/Mosaic_Years_Since_Greatest_Change_2010.img")
#   Change_rate <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_Change_rate_2010.img")
#   Change_persistence <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_Change_persistence.img")
#   Change_magnitude <- raster("G:/ForGiani_WorkingImputation/data/change/Mosaic_Change_magnitude_variation_2010.img")
#   ##load topographic metrics
#   Elevation <- raster("G:/ForGiani_WorkingImputation/data/topo/CDED.img")
#   Slope <- raster("G:/ForGiani_WorkingImputation/data/topo/SLOPE_DEG.img")
#   TWI <- raster("G:/ForGiani_WorkingImputation/data/topo/TWI.img")
#   TSRI <- raster("G:/ForGiani_WorkingImputation/data/topo/TSRI.img")
  
  TC <- stack(file.path(TC_dir, paste(zone,"_TC_Mosaic_2010.img",sep = ''),fsep = .Platform$file.sep))
  PreChange_persistence <- raster(file.path(change_dir, "Mosaic_PreChange_persistence_2010.img", fsep = .Platform$file.sep))
  PreChange_magnitude <- raster(file.path(change_dir, "Mosaic_PreChange_magnitude_variation_2010.img", fsep = .Platform$file.sep))
  PreChange_evolution <- raster(file.path(change_dir, "Mosaic_PreChange_evolution_rate_2010.img", fsep = .Platform$file.sep))
  PostChange_persistence <- raster(file.path(change_dir, "Mosaic_PostChange_persistence_2010.img", fsep = .Platform$file.sep))
  PostChange_magnitude <- raster(file.path(change_dir, "Mosaic_PostChange_magnitude_variation_2010.img", fsep = .Platform$file.sep))
  PostChange_evolution <- raster(file.path(change_dir, "Mosaic_PostChange_evolution_rate_2010.img", fsep = .Platform$file.sep))
  Years_Since_Greatest_Change <- raster(file.path(change_dir, "Mosaic_Years_Since_Greatest_Change_2010.img", fsep = .Platform$file.sep))
  Change_rate <- raster(file.path(change_dir, "Mosaic_Change_rate_2010.img", fsep = .Platform$file.sep))
  Change_persistence <- raster(file.path(change_dir, "Mosaic_Change_persistence.img", fsep = .Platform$file.sep))
  Change_magnitude <- raster(file.path(change_dir, "Mosaic_Change_magnitude_variation_2010.img", fsep = .Platform$file.sep))
  Elevation <- raster(file.path(topo_dir, "CDED.img", fsep = .Platform$file.sep))
  Slope <- raster(file.path(topo_dir, "SLOPE_DEG.img", fsep = .Platform$file.sep))
  TWI <- raster(file.path(topo_dir, "TWI.img", fsep = .Platform$file.sep))
  TSRI <- raster(file.path(topo_dir, "TSRI.img", fsep = .Platform$file.sep))
  
  
  
  beginCluster(30,type="SOCK")
  ##make a rasterstack of all explantory variables
  env.vars.stack <- stack(TC,
                          PreChange_persistence, PreChange_magnitude, PreChange_evolution,
                          PostChange_persistence, PostChange_magnitude, PostChange_evolution,
                          Years_Since_Greatest_Change, Change_rate, Change_persistence, Change_magnitude,
                          Elevation, Slope, TWI, TSRI)
  
  ## rename layers in raster stack to match variable names in RF model
  names(env.vars.stack) <- c("TCB","TCG","TCW","TCA","TCD",
                             "PreChange_persistence","PreChange_magnitude","PreChange_evolution",
                             "PostChange_persistence","PostChange_magnitude","PostChange_evolution",
                             "Years_Since_Greatest_Change", "Change_rate","Change_persistence","Change_magnitude",
                             "Elevation","Slope","TWI","TSRI")
  
  
  
  ##Crop and mask raster stack to just the bnd area of interest and forest class.
  ##saves processing time if your raster data is much larger than your area of interest.
  
  env.vars.stack.2 <- crop(env.vars.stack, bnd, progress = "text")  # SUPER SLOW !!!
  env.vars.stack.2 <- mask(env.vars.stack.2, fc.mask, progress = "text")  # sk.fc.mask.2 NON EXISTENT
  
  
  ##wrapper function for previously created yai model either gnn or random forest
  ##converts the IDs from character to numeric.
  predfun.yai<-function(yaimodel, data){
    rownames(data)<-paste("X",1:nrow(data),sep = "") ## so that newtargets doesn't get confused by duplicate row IDs.
    return(as.numeric(newtargets(yaimodel,data,k = 1)$neiIdsTrgs))
  }
  
  
  # load("H:/SK_nn_forest/sk_models/rf_gnn_models/yai_rf_gnn_sk.RData")  # ALREADY LOADED ABOVE
  # load("H:/SK_nn_forest/Geordie_wkg/wkg/imputation results/en_vars_stack2.R")
  ##map it using clusterR
  ##also spit processing time
  
  # showConnections(all=T)
  # closeAllConnections()
  
  start.time <- Sys.time()
  # ##beginCluster()
  # detectCores()
  # ## [1] 4
  # # Create cluster with desired number of cores
  # cl <- makeCluster(24)
  # # Register cluster
  # registerDoParallel(cl)
  # # Find out how many cores are being used
  # getDoParWorkers()
  
  # CHECK VARIABLES' NAME MATCH
  # prediction on 883'882'701 samples (883 millions!!!), all over Canada will be: 400'000'000 ha / 0.09 ha/pix = 4'444'444'444 pixels (4 billion pixels)
  # equivalent to calling: predict(env.vars.stack.2, yai.rf.sk.500, fun=predfun.yai) with predfun.yai which in turn calls newtargets() as main predicting function
  # "filename", "progress", "ext", "overwrite" and "m" are all arguments of the predict function
  clusterR(env.vars.stack.2, predict, args=list(yai.rf.sk.500, fun=predfun.yai), filename=file.path(results_dir, "sk_rf_2010_fcid.img", fsep = .Platform$file.sep), progress="text", ext=extent(bnd), overwrite=T, m=2000)
  
  endCluster()
  
  # predict(env.vars.stack.2, yai.rf.sk.500, filename = "h:/SK_nn_forest/Geordie_wkg/wkg/imputation results/sk_subset_rf_2010_fcid.img",progress = "text", ext = extent(bnd), overwrite = T, m = 2000)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  write.csv(time.taken, file = "time.taken.csv")
  warnings()
  
  ##load imputation raster of fcids for sk_gnn_2010
  sk.rf.2010.fcid <- raster(file.path(results_dir, "sk_rf_2010_fcid.img", fsep = .Platform$file.sep))
  
  ##create ID column from rownames
  Y.trn.attach$ID <- as.numeric(rownames(Y.trn.attach))
  
  ##ratify the raster and set up a lookup table
  r <- ratify(sk.rf.2010.fcid)
  rat <- levels(r)[[1]]
  
  ##create a lookup table from the ancillary dataframe
  lookup <- merge(rat, Y.trn.attach, by.x="ID", by.y = "ID")
  
  ##attach lookup table to raster
  levels(r) <- lookup
  
  ##deratify to create a single layer biomass map
  biomass <- deratify(r, 'total_biomass') 
  dataType(biomass)
  writeRaster(biomass, file.path(results_dir, "sk_rf_2010_total_biomass.img", fsep = .Platform$file.sep), progress="text",format='HFA', dataType= 'FLT4U',overwrite=TRUE) 
  
  ##delete temp folder with temp rasters
  unlink(raster_tmp_dir, recursive = T, force = T)
  
  # clock UTM zone time
  temp.toc <- proc.time()-temp.tic[3]
  print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))
  
}

# stopCluster(cl)

# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog4, total elapsed time:",seconds_to_period(toc[3])))
