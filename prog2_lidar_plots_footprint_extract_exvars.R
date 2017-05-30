## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog2_lidar_plots_footprint_extract_exvars.R                             
## Objective: Extract values for all explanatory rasters spatially coincident with sampled lidar plot polygons for both training and validation plots
##            Mean plot level values of explanatory variables is weighted average based on relative proportion of each pixel within the polygon (normalized weights)

#### TO DO -------------------------------------------------------------------

## STILL TO DO:

# Prior to actual run:
# - load correct param_file
# - set params2$copy.data <- T
# - check params2$fitted.data
# - check copying on .dat and not on .hdr check ## TO TEST SCRIPT in Functions_NatMapping_R.R
# - uncomment unlink() to delete temp
# - use foreach on blocks
# - use params2$contin.sampling.shp <- 'poly' for sampling
# - check no redefinition of paramsGL$zones 
# - check parameters
# - delete all UTMzones folders in Landsat_dir E:\NTEMS otherwise copy_paste_UTMdata() will not be run
# - uncomment final unlink() to delete copied folders

## SOLVED:
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
# -V no point in saving a dataset (table) for mapping, do a script creating variables on the fly for prediction/mapping -- yes, thats the way to do it, consider saving a TC components file (5 layers) with the new narrow extent
# -V add ArcGIS/ArcPy AreaSolarRadiation() ? -- no, super-slow to compute and already included in TSRI
# -v check extents (Topo by Geordie vs rest) -- Geordie has a script to change and resave all NTEMS layers at the new narrow extent (less overlap) with the same folder structure
# -V convert "Year of " variables to "Years since" to avoid weird averages
# -V adapt behavior of params2$pred.names.lg and params2$pred.names.sh wrt list vs vector
# -V how to ignore unclassified in Change_attribution? -- we'll use it with all the classes separated at the moment
# -V Keep following predictors? 
#     disturbed -- No, info included in "First_Change_Year"
#     Years_Since_Greatest_Change, Years_Since_Last_Change  Years_Since_First_Change -- Yes, included instead of "_Change_Year" variables to be able to take averages
#     Aspect, Aspect_cos -- no, included in TSRI, TWI
# -V Differences in the rasters read here between Harold's files (eg. "Mosaic_PreChange_persistence_2010.img") and NTEMS ones 
#   (eg."SRef_13S_PreChange_persistence.dat"). Apparently not only in terms of the extent but also in terms of pixel values. -- Not to pay attention to bc layers have been recomputed since when Harlod did his study
# -V check metric for disturbed and undisturbed polygons... to do (all these histograms that are now removed)? -- check instead histograms for all UTM zone separately, scatterplots, etc.
# -V climate variables at 1km (finer res?) -- Resolution is 1KM, Doug is not sure precipitation will help, Temp probably will (for prediciting high trees) and he's investigating how to integrate these layers (12 monthly measurements over several years), wait for his results.
# -V feasible to recompute on the fly the new variables for mapping phase? -- Yes, but let's try to stick to raw variables as much as possible
# -V yr of to yrs since: params2$min.year.of <- 1984, params2$max.year.of <- 2011 ?? -- No, better to keep Year_of, sample with NN and plot relationships only for year_of != 0
# -V use nearest neighbor for sampling years layers (no average) -- OK like that
# -V TSRI different values the two different layers -- layers have been recomputed since Harold's study
# -V hours of daylight from lat long -- only a function of Lat, annual averages out, length of day at Aug 1st is better but ultimately is not used bc highly correlated with latitude
# -V #### TOPO KO  -- split extract for each type of predictor, so now Topo OK
# -V check new Ch_attr -- now recode is OK (merge all unclassified)
# -V NTEMS cut to different extents -- separate extract() for Topo
# -V rastertemp directory is used? -- dont know but doesnt hurt to have it there
# -V check projections to ensure we're sampling proper locations if-else-stop blocks inserted in script
# -V uncomment params2$pred.names.sh -- now extract for topo is separate so topo is included in pred.names 
# -V test with new ceiling() in parallel extract on blocks of poly.val -- added if-else-stop block that checks for equal nr of rows of the input and outputs of the foreach

#### READS/WRITES ------------------------------------------------------------

## READS:
# - "<UTMzone>_pt_centerpt_training_validation.shp": (from prog1) point shp with centerpoints of selected polygons
# - "<UTMzone>_poly_training_validation.shp": (from prog1) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)

## WRITES:
# - "poly_training_validation_exvars_extract.csv": (for prog3a and prog3) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

#### INIT --------------------------------------------------------------------

print('Prog2: explanatory variables extraction')

rm(list=ls()) ## clear all variables

# param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg_BOREAL/AllUTMzones_paramsGL.Rdata"
param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")

#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params2 <- list()

# params2$copy.data <- F   ## wheter or not to copy the data from remot disks (set to True for actual run)
params2$copy.data <- T

# XXXXXXXXXXX TODEL XXXXXXXXXX
# paramsGL$zones <- c("UTM10S", "UTM12S")
# XXXXXXXXXXXXXXXXXXXXXXXXXXXX

## choice of sampling unit for continuous variables (bands and topo)
params2$contin.sampling.shp <- 'poly'
# params2$contin.sampling.shp <- 'cpt'

# params2$fitted.data <- F
params2$fitted.data <- T

params2$fitted.proxies.dir <- file.path(Landsat_dir, 'UTM_results', '2017_02_27_SAMPLES_treed_changes_plots', fsep = .Platform$file.sep)

params2$lidar.sources <- c("BOREAL", "NONBOREAL")
# params2$lidar.sources <- c("NONBOREAL")

params2$transect.names <- c("AFRF", "CR04", "CR08", "MK10", "Quesnel", "Tofino")
params2$transect.years <- c(2008, 2004, 2008, 2010, 2008, 2005)
# params2$transect.names <- c("CR04", "CR08")
# params2$transect.years <- c(2004, 2008)

## 22 predictors (.lg names are generally used to load files, .sh are shorter dataframe column names):

## 6 input spectral bands
bands.names.lg <- list('band1', 'band2', 'band3', 'band4', 'band5', 'band7')
bands.names.sh <- list('b1', 'b2', 'b3', 'b4', 'b5', 'b7')

## 5 derived Tasseled Cap components (added in the very end)
TC.names.lg <- list('TC_Brightness', 'TC_Greenness', 'TC_Wetness', 'TC_Angle', 'TC_Distance')  
TC.names.sh <- list('TCB', 'TCG', 'TCW', 'TCA', 'TCD')

## 3 Vegetation indices
VI.names.lg <- list('Normalized Difference Vegetation Index', 'Enhanced Vegetation Index', 'Normalized Burn Ratio')  
VI.names.sh <- list('NDVI', 'EVI', 'NBR')

## 1 unaverageable "years of/since" change variables (separate extract())
year.cng.input.names.lg <- list("Greatest_Change_Year")
year.cng.names.lg <- list("Years_Since_Greatest_Change")
year.cng.names.sh <- list("YrsSince_GrCh")

## 1 unaverageable categorical layer of type of changes
cngattr.names.lg <- list("Change_attribution")
cngattr.names.sh <- list("Ch_attr")

## 4 topography related variables
topo.names.lg <- list("Elevation", "Slope", "Topographic_wetness_index", "Topographic_solar_radiation_index")
topo.names.sh <- list("Elev", "Slope", "TWI", "TSRI")

## 2 spatially contiuous variables for de-trending: coordinates
coords.names.lg <- list("Longitude", "Latitude")
coords.names.sh <- list("Long", "Lat")

params2$pred.names.lg <- list(bands=bands.names.lg, TC=TC.names.lg, VI=VI.names.lg, yrs=year.cng.names.lg, cngattr=cngattr.names.lg, topo=topo.names.lg, coords=coords.names.lg)
params2$pred.names.sh <- list(bands=bands.names.sh, TC=TC.names.sh, VI=VI.names.sh, yrs=year.cng.names.sh, cngattr=cngattr.names.sh, topo=topo.names.sh, coords=coords.names.sh)

params2$years.since.no.change <- 50  ## value to assign to pixels with 0 as Year_of_GreatCh (no change areas are set as last changed in this many years ago)

params2$disk.names <- list("//frst-frm-2234b/I", "//frst-frm-2234b/K")

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep)
save(params2, file = param_file_prog2)

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "dplyr",    ## to be loaded before foreach to avoid "assertion failed" errors
                      "landsat",
                      "lattice",
                      "lubridate", 
                      "doParallel", 
                      "foreach", 
                      "data.table"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=T)
}



#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time

temp.dir <- file.path(Landsat_dir, "temp", fsep = .Platform$file.sep)

idx.log <- 1
copy.time.UTMzone <- vector("list", length(paramsGL$zones)+length(params2$transect.names)) 
full.df.CAN.list <- vector("list", length(params2$lidar.sources)) 
for (ls in 1:length(params2$lidar.sources)) {
  
  lidar.source <- params2$lidar.sources[ls]
  
  if (lidar.source == "BOREAL") {
    zones <- paramsGL$zones
  } else {
    zones <- params2$transect.names
  }

  full.table <- vector("list", length(zones)) 
  for (z in 1:length(zones)) {
    
    temp.tic <- proc.time() ## start clocking time for each UTM zone
    
    zone <- zones[z]
    if (lidar.source == "BOREAL") {
      zone.nr <- substr(zone, 4, nchar(zone))
      year <- 2010
    } else {
      zone.nr <- "10S"
      year <- params2$transect.years[z]
    }
    print(paste('Processing zone', zone))   ## converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
    wkg_dir = file.path(base_wkg_dir, zone, fsep = .Platform$file.sep)
  
  #### COPY DATA FROM REMOTE DISKS ---------------------------------------------
    
    ## call copy_paste_UTMdata() only if we specify to copy data
    if (params2$copy.data) {
      
      ## first delete folder after having read the data and make room for the next zone's data
      unlink(temp.dir, recursive = T, force = T)
      
      log.copy.paste <- copy_paste_UTMdata(params2$disk.names, zone, zone.nr, year, year.cng.input.names.lg, temp.dir)
      if (length(log.copy.paste$list.missing.files)>0) {
        stop(sprintf("Missing files for zone %s", zone))
      }
      copy.time.UTMzone[idx.log] <- paste(zone, "copy time:", seconds_to_period(log.copy.paste$toc[3]))
      idx.log <- idx.log+1
    }
    
    ## running large rasters results in the hard drive filling up with temp gribs
    ## deal with it by telling raster package to write temp files in a new defined directory within the working directory
    ## then remove the directory at the end of session
    raster_tmp_dir <- "raster_tmp"
    dir.create(raster_tmp_dir, showWarnings = F, recursive = T)  ## creates temp directory in BAP_Imputation_working/wkg/<UTMzone>
    rasterOptions(tmpdir = raster_tmp_dir)
  
  #### READ FILES -----------------------------------------------------
    
    ## read training and validation polygons and centerpoints
    poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone, "_poly_training_validation", sep = ''))
    cpt.poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone, "_pt_centerpt_training_validation", sep = ''))
  
    ## Landsat bands for the year of interest (6 layers), TC layer will be computed later (first for TRN/VAL and later for all Canada when mapping)
    if (!params2$fitted.data) {
      ## Regular rasters
      cmd <- sprintf('Bands <- brick(file.path("%s", "UTM_%s", "Results", "proxy_values", "SRef_UTM%s_%d_proxy_v2.dat", fsep = .Platform$file.sep))', temp.dir, zone.nr, zone.nr, year)  
    } else {
      if (lidar.source == "BOREAL" ) {
        ## Fitted 2010 data for BOREAL transect
        cmd <- sprintf('Bands <- brick(file.path("%s", "plots_crown", "proxies", "UTM%s_proxy_plots_crown.dat", fsep = .Platform$file.sep))', params2$fitted.proxies.dir, zone.nr)  
      } else {
        if (zone == 'CR08') {
          ## Fitted 2008 data for CR dataset
          cmd <- sprintf('Bands <- brick(file.path("%s", "plots_crown_10", "proxies", "UTM10S_proxy_plots_crown_CR08.dat", fsep = .Platform$file.sep))',  params2$fitted.proxies.dir)  
        } else {
          ## Fitted data for all the other NONBOREAL datasets (other years)
          cmd <- sprintf('Bands <- brick(file.path("%s", "plots_crown_10", "proxies", "UTM10S_proxy_plots_crown_noCR08.dat", fsep = .Platform$file.sep))',  params2$fitted.proxies.dir)  
        }
      }
    }
    eval(parse(text=cmd))
    
    ## "Greatest_Change_Year" (1 layer)
    ## (formerly read all change raster layers from folder <UTM_zone.nr>\Results\Change_metrics\
    ## temp.dir goes at the end as a variable bc if left inside sprintf as a string is not detected by foreach and thus not passed to the cores (error: "object 'temp.dir' not found")
    for (var in year.cng.input.names.lg) {
      cmd <- sprintf('%s <- raster(file.path("%s", "UTM_%s", "Results", "Change_metrics", "SRef_%s_%s.dat", fsep = .Platform$file.sep))', var, temp.dir, zone.nr, zone.nr, var)  
      eval(parse(text=cmd))
    }
    
    ## change attribution (1 layer)
    cmd <- sprintf('Change_attribution <- raster(file.path("%s", "UTM_%s", "Results", "change_attribution", "UTM%s_Change_attribution_complete.dat", fsep = .Platform$file.sep))', temp.dir, zone.nr, zone.nr)  
    eval(parse(text=cmd))
    
    ## topographic metrics (4 layers)
    Elevation <- raster(file.path(topo_dir, "DEM", paste("UTM", zone.nr, "_DEM.dat", sep= ''), fsep = .Platform$file.sep))
    Slope <- raster(file.path(topo_dir, "TopoMetrics", paste("UTM", zone.nr, "_TopoMetrics.dat", sep= ''), fsep = .Platform$file.sep), band = 1)  ## band 1 is Slope whereas band 2 and 3 are Aspect and Shaded Relief
    Topographic_wetness_index <- raster(file.path(topo_dir, "TWI", paste("UTM", zone.nr, "_TWI.dat", sep= ''), fsep = .Platform$file.sep))
    Topographic_solar_radiation_index <- raster(file.path(topo_dir, "TSRI",  paste("UTM", zone.nr, "_TSRI.dat", sep= ''), fsep = .Platform$file.sep))
  
  #### CREATE RASTER STACKS ---------------------------------------------------
  
    ## separate read of continuous variables (bands and topo) since we don't have topo and other NTEMS data matching in extent
    vars.bands <- c(bands.names.sh)   ## to rename columns of final dataframe
    str.bands <- paste(c('Bands'), collapse=", ")    ## creates string of variable names to input in stack()
    vars.topo <- c(topo.names.sh)  
    str.topo <- paste(c(topo.names.lg), collapse=", ")
    vars.years <- year.cng.names.sh
    str.years <- paste(year.cng.input.names.lg, collapse=", ")
    vars.categ <- cngattr.names.sh
    str.categ <- paste(cngattr.names.lg, collapse=", ")
    
    cmd <- sprintf("exvars.stack.bands <- stack(%s)", str.bands)  ## stack() already checks if projections of all the layers match
    eval(parse(text=cmd))
    names(exvars.stack.bands) <- vars.bands
    
    cmd <- sprintf("exvars.stack.topo <- stack(%s)", str.topo)  
    eval(parse(text=cmd))
    names(exvars.stack.topo) <- vars.topo
    
    cmd <- sprintf("exvars.stack.years <- stack(%s)", str.years)  
    eval(parse(text=cmd))
    names(exvars.stack.years) <- vars.years  ## already change here the name from "Year of" to "Years since"
    
    cmd <- sprintf("exvars.stack.categ <- stack(%s)", str.categ) 
    eval(parse(text=cmd))
    names(exvars.stack.categ) <- vars.categ
    
    ## check if coordinate systems are the same between raster stack and shapefiles
    if ( !all(sapply(list(proj4string(poly.training.validation), proj4string(cpt.poly.training.validation), proj4string(exvars.stack.topo), proj4string(exvars.stack.years), proj4string(exvars.stack.categ)), FUN = identical, proj4string(exvars.stack.bands))) ) {
      stop(sprintf("%s: projections of raster stack and shapefiles loaded in prog2 do not match.", zone))
    }
  
  #### PARALLEL LOOP ON BLOCKS OF poly.training.validation ----------------------
  
    sizeBlocks <- ceiling( nrow(poly.training.validation@data) / min(detectCores(), nrow(poly.training.validation@data)) )  ## to avoid problems when number of rows of the matrix is smaller than the nr of cores 
    nblocks <- ceiling(nrow(poly.training.validation@data)/sizeBlocks)   ## to avoid problems when sizeBlocks is small (nblocks will allow to cover just a little more, or exactly, the size of data)
    nr.clusters <- nblocks
    cl <- makeCluster(nr.clusters)
    registerDoParallel(cl)
    exvars.extract <- foreach(blck = 1:nblocks, .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=T for more info when debugging
      ## WHEN NOT USING FOREACH, TESTING PHASE
    # nblocks <- 19
    # for (blck in 17:nblocks) {
    #   sizeBlocks <- 1
      ## WHEN NOT USING FOREACH, TESTING PHASE
      
      indPart <- seq(from=(blck-1)*sizeBlocks+1, to=min(blck*sizeBlocks, nrow(poly.training.validation@data)), by=1)
  
  #### EXTRACT VALUES OVER POLYGONS ---------------------------------------------------
      
      ## selection of the sampling unit for continuous variables
      if (params2$contin.sampling.shp == 'poly') {   
        sampling.shp <- poly.training.validation[indPart, ]   ## for actual study we use polygons
      } else if (params2$contin.sampling.shp == 'cpt') { 
        sampling.shp <- cpt.poly.training.validation[indPart, ]    ## to speed up tests read only centerpoint
      }
      
      ## continuous variables: use 'poly' in the option above for actual run
      ## extract mean of exvars weighted by nomralized proportion contribution of each cell in polygon (weights within each polygon sum to 1)
      ## automatically assigns an ID starting at 1 block of samples assigned to the cores (to be removed later on bc meaningless)
      exvars.extract.bands <- extract(exvars.stack.bands, sampling.shp, fun = mean, na.rm = F, weights = T, normalizeWeights = T,
                            cellnumbers = F, small = F, df = T, factors = F, sp = F)  # extract is a function that actually accesses the values of the rasters, MEMORY ISSUES?
      exvars.extract.topo <- extract(exvars.stack.topo, sampling.shp, fun = mean, na.rm = F, weights = T, normalizeWeights = T,
                                     cellnumbers = F, small = F, df = T, factors = F, sp = F) 
  
      ## categorical/un-averageable variables: use cpt.poly.training.validation so that we extract only one single value, that of the pixels right below the center of the plot
      exvars.extract.years <- extract(exvars.stack.years, cpt.poly.training.validation[indPart, ], na.rm = F,
                         cellnumbers = F, small = F, df = T, factors = F, sp = F)
      exvars.extract.years$ID <- NULL
      
      ## to transform from "Year of" to "Years since"
      idx.zeros <- exvars.extract.years$YrsSince_GrCh == 0    ## indices of pixels that never changed (here exvars.extract.years$YrsSince_GrCh is still the "Year of Change", and 0 means never changed)
      exvars.extract.years$YrsSince_GrCh[idx.zeros] <- params2$years.since.no.change  ## start converting to YrsSince_GrCh: set large value (50) for pixels that never changed
      exvars.extract.years$YrsSince_GrCh[!idx.zeros] <- year - exvars.extract.years$YrsSince_GrCh[!idx.zeros]  ## for the rest of pixels, set YrsSince_GrCh as the year of acquisition (year) minus the "Year of Change" 
      idx.not.yet.changed <- exvars.extract.years$YrsSince_GrCh < 0   ## retrieve indices of pixels that have changed after the LiDAR acquisition year (will now have a negative value because "Year of Change" was larger than year)
      exvars.extract.years$YrsSince_GrCh[idx.not.yet.changed] <- params2$years.since.no.change  ## for pixels having changed after our TARGET_YEAR (2010) or the year of acquisition (for BC data) let's assign the no-change value (50)
      
      exvars.extract.categ <- extract(exvars.stack.categ, cpt.poly.training.validation[indPart, ], na.rm = F,
              cellnumbers = F, small = F, df = T, factors = F, sp = F)
      exvars.extract.categ$ID <- NULL
      exvars.extract.categ$Ch_attr[idx.not.yet.changed] <- 0   ## set type of change to "No change" for not yet changed pixels
      
      cbind(exvars.extract.bands, exvars.extract.topo, exvars.extract.years, exvars.extract.categ)
  
      ## WHEN NOT USING FOREACH, TESTING PHASE
      # exvars.extract <- cbind(exvars.extract.bands, exvars.extract.topo, exvars.extract.years, exvars.extract.categ)
      ## WHEN NOT USING FOREACH, TESTING PHASE
      
    }  ## end for on nblocks
    
    stopCluster(cl)
    
    ## check if nr of rows is the same in input and output of the foreach (to avoid issues with block indexes)
    if ( nrow(poly.training.validation@data) != nrow(exvars.extract) ) {
      stop(sprintf("Zone %s: nr. of rows is different for poly.training.validation and exvars.extract", zone))
    }
    
    ## delete temp folder with temp rasters
    unlink(raster_tmp_dir, recursive = T, force = T)
  
  #### ADD LAT/LONG & EXTRA INFO TO DATA FRAME IN LOOP  ----------------------------------
  
    ## coordinates of plot polygons
    long.lat <- as.data.frame(coordinates(spTransform(cpt.poly.training.validation, CRS("+proj=longlat"))))
    coords <- as.data.frame(long.lat)
    colnames(coords) <- coords.names.sh
    
    UTMzone.ID <- as.data.frame(rep(zone.nr, nrow(cpt.poly.training.validation)))
    colnames(UTMzone.ID) <- "UTMzone"
    
    ## bind training polygon information with UTMzone ID, extracted explanatory variables and coordinates
    infos <- cpt.poly.training.validation@data[,c("FCID","POLY250ID","TV")]
    full.table[[z]] <- cbind(infos, UTMzone.ID, exvars.extract, coords)
    
    ## clock UTM zone time
    temp.toc <- proc.time()-temp.tic[3]
    print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))
    
  }  ## end for on length(zones)
  
  full.df.CAN.list[[ls]] <- as.data.frame(rbindlist(full.table))
  
}  ## end for on params2$lidar.sources

full.df.CAN <- as.data.frame(rbindlist(full.df.CAN.list))

#### ADD SPECTRAL INDICES, RECODE CH_ATTR OUTSIDE LOOP & WRITE CSV ------------------------------------

all.bands <- as.matrix(full.df.CAN[, c("b1", "b2", "b3", "b4", "b5", "b7")]/10000)  ## to bring back the values to reflectance (needed for EVI)

## compute TCs
TC.list <- bands_2_TC(all.bands)
TC.df <- cbind(TC.list$TCs, TC.list$TCang, TC.list$TCdist)
colnames(TC.df) <- TC.names.sh

## compute VI
NDVI <- (all.bands[,"b4"]-all.bands[,"b3"])/(all.bands[,"b4"]+all.bands[,"b3"])
EVI <- 2.5*( (all.bands[,"b4"]-all.bands[,"b3"]) / (1 + all.bands[,"b4"] + 6*all.bands[,"b3"] - 7.5*all.bands[,"b1"]) )
NBR <- (all.bands[,"b4"]-all.bands[,"b7"])/(all.bands[,"b4"]+all.bands[,"b7"])
VI.df <- data.frame(NDVI, EVI, NBR)
colnames(VI.df) <- VI.names.sh

## recode Change attr from 10 to 5 classes:
## original 10 classes:  0="No change", 1="Fire", 2="Harvesting", 3="Non-stand replacing", 4="Road", 5="Unclassified (Fire)", 6="Unclassified (Harvesting)", 7="Unclassified (Non-stand replacing)", 8="Unclassified (Road)", 9="Unclassified (Well-site)"
## retained 5 classes:   0="No change", 1="Fire", 2="Harvesting", 3="Non-stand replacing", 4="Infrastructure"
Ch_attr.10classes <- as.matrix(full.df.CAN[, "Ch_attr"])  
Ch_attr.reclassified <- Ch_attr.10classes
Ch_attr.reclassified[Ch_attr.10classes == 5] <- 1  ## merge Fire and unclassified Fire
Ch_attr.reclassified[Ch_attr.10classes == 6] <- 2  ## merge Harvesting and unclassified Harvesting
Ch_attr.reclassified[Ch_attr.10classes == 7] <- 3  ## merge Non-stand replacing and unclassified Non-stand replacing
Ch_attr.reclassified[Ch_attr.10classes == 8] <- 4  ## merge Road and unclassified Road
Ch_attr.reclassified[Ch_attr.10classes == 9] <- 4  ## merge Road and unclassified Well-site, now called Infrastructure

## substitute the old Ch_attr column with the newly recoded one
full.df.CAN[, "Ch_attr"] <- as.data.frame(Ch_attr.reclassified)

## add the 5 new TC components and the 3 VIs
full.df.CAN <- cbind(full.df.CAN, TC.df, VI.df)

## reorder columns to match predefined order and remove ID columns
full.df.CAN <- full.df.CAN[, c("FCID", "POLY250ID", "TV", "UTMzone", unlist(params2$pred.names.sh))]

## write full.df.CAN (former "poly.training.validation.exvars.extract") as csv after sortrows (arrange) based on FCID (for direct comparison with lidar_metrics_mean_training_validation.csv)
full.df.CAN$FCID <- as.character(full.df.CAN$FCID)
if (!params2$fitted.data) {
  write.csv(arrange(full.df.CAN, FCID), file = file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", sep = ''))
} else {
  write.csv(arrange(full.df.CAN, FCID), file = file.path(base_wkg_dir, "poly_training_validation_exvars_extract_fitted.csv", sep = ''))
}

## ultimately delete the temp folder
unlink(temp.dir, recursive = T, force = T)


#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog2, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)

## write log file
unlink(file.path(base_wkg_dir, "log_prog2.txt", sep = ''), recursive = T, force = T)
list.to.print <- append(copy.time.UTMzone, end.message)
lapply(list.to.print, write, file.path(base_wkg_dir, "log_prog2.txt", sep = ''), append=T)

