#### CODE INFOS -------------------------------------------------------------------

## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog4_check_maps.R                            
## Objective: Check that predicted values on maps match those predicted on plots

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - use foreach on blocks
# - check no redefinition of paramsGL$zones and subsetting of poly.validation: look for ##---
# - check parameters

# - check script for new boxplots v3
# - check influence of rasterOptions(maxmemory = 1e+09), rasterOptions(chunksize = 1e+08)


## SOLVED:
# -V change back and remove _v2 -- now we have separate folders for each version of the csvs
# -V save separate correlation trends plots -- plots are now produced


#### READS/WRITES ------------------------------------------------------------

## READS:
# - "<UTMzone>_pt_centerpt_training_validation.shp": (from prog1) point shp with centerpoints of selected polygons
# - "<UTMzone>_poly_training_validation.shp": (from prog1) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)

## WRITES:
# - "poly_training_validation_exvars_extract.csv": (for prog3a and prog3) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

#### INIT --------------------------------------------------------------------

print('Prog4: check maps')

rm(list=ls()) ## clear all variables

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params2.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params3.Rdata"
load(param_file)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")

## Update results directory to point to the FINAL results one
base_results_dir <- file.path(base_dir, "FINAL_RUN", "results", fsep=.Platform$file.sep)

#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

##-------
paramsGL$zones <- c("UTM12S", "UTM13S")
# paramsGL$zones <- c("UTM13S")

##-------

params4 <- list()

params4$subsetting <- F    ## to subset the full map csvs (for debugging in development phase)
params4$check.plot.vs.pix <- T   ## wheter to run the plot vs pixel values comparison, section 1)
params4$check.temporal.val.traj <- F  ## wheter to run the single pixel trajectory analysis on val set, section 2)
params4$check.map.coherence <- F  ## wheter to run the map coherence analysis, section 5)
params4$mapped.years <- seq(from=1984, to=2012)
# params4$temp.sampling.shp <- 'poly'
params4$temp.sampling.shp <- 'cpt'
# params4$incoher.quant <- c(0.3, 0.7)
params4$lims.el_p95 <- c(10, 12)
params4$lims.p_1r_2m <- c(40, 60)
params4$targ.to.plot <- "el_p95"
params4$nr.sampled.pix <- 500
params4$ecozone.codes <- c(4, 6, 18)
params4$ecozone.names <- c("Boreal Plains", "Boreal Shield West", "Taiga Shield West")
# params4$disk.names <- list("//frst-cdw-2231j/I", "//frst-cdw-2231j/J")

params4$multipl <- list(elev_mean=1000, elev_stddev=1000, elev_cv=1000, elev_p95=1000, 
                        percentage_first_returns_above_2m=100, percentage_first_returns_above_mean=100,
                        loreys_height=1000, basal_area=100, gross_stem_volume=100, total_biomass=10)

param_file_prog4 = file.path(base_wkg_dir, 'AllUTMzones_params4.Rdata', fsep = .Platform$file.sep) 
save(params4, file = param_file_prog4)

#### LOAD PACKAGES ----------------------------------------------------------

## install data.table version with fwrite
# install.packages("data.table", repos = "https://Rdatatable.github.io/data.table", type = "source")

list.of.packages <- c("ggplot2",
                      "gridExtra",
                      "reshape2",
                      "rgdal",
                      "raster",
                      "GSIF",
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

## running large rasters results in the hard drive filling up with temp gribs
## deal with it by telling raster package to write temp files in a new defined directory within the working directory
## then remove the directory at the end of session (default is in "C:\Users\gmatasci\AppData\Local\Temp\RtmpUFQ26C\raster")
raster_tmp_dir <- file.path(base_wkg_dir, "raster_tmp", fsep = .Platform$file.sep)
dir.create(raster_tmp_dir, showWarnings = F, recursive = T)
rasterOptions(tmpdir=raster_tmp_dir)

## If we check for plot/pixel correspondence load predictions on Val set 
if (params4$check.plot.vs.pix | params4$check.temporal.val.traj) {  
  plot.based.predictions <- readOGR(dsn = base_results_dir, layer = 'val_predictions_RF')
  colnames(plot.based.predictions@data)[ncol(plot.based.predictions@data)] <- "FCID"
  plot.based.predictions@data$FCID <- as.character(plot.based.predictions@data$FCID) 
}

full.table.targets <- vector("list", length(paramsGL$zones)) 
full.table.years <- vector("list", length(paramsGL$zones)) 
for (z in 1:length(paramsGL$zones)) {
  
  temp.tic <- proc.time() ## start clocking time for each UTM zone
  
  zone <- paramsGL$zones[z]
  zone.nr <- substr(zone, 4, nchar(zone))
  print(paste('Processing zone', zone))   ## converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)

#### READ FILES AND CREATE RASTER STACKS ----------------------------------------
  
  ## read training and validation shapefiles: polygons and centerpoints
  if (params4$check.plot.vs.pix | params4$check.temporal.val.traj) { 
    poly.training.validation <- readOGR(dsn = wkg_dir, layer = paste(zone,"_poly_training_validation",sep = ''))
    cpt.poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_pt_centerpt_training_validation",sep = ''))
  }
  
  ## read 2010 maps for all targets
  if (params4$check.plot.vs.pix | params4$check.map.coherence) {
    
    for (targ in params3$targ.names.lg) {
      cmd <- sprintf('%s <- raster(file.path("%s", "UTM_results", "UTM_%s", "RF", "UTM_%s_%s_RF_2010.dat", fsep = .Platform$file.sep))', targ, Landsat_dir, zone.nr, zone.nr, targ)  
      eval(parse(text=cmd))
    }
    
    vars.targets <- c(params3$targ.names.sh)  
    str.targets <- paste(c(params3$targ.names.lg), collapse=", ")
    
    cmd <- sprintf("stack.targets <- stack(%s)", str.targets)  
    eval(parse(text=cmd))
    names(stack.targets) <- vars.targets
    
  }
  

#### EXTRACT ALL RESP. VARIABLES VALUES FROM 2010 MAPS ----------------
  
  if (params4$check.plot.vs.pix) {
    
    ## check if coordinate systems are the same between raster stack and shapefiles
    if ( !all( sapply( list(proj4string(poly.training.validation)), FUN = identical, proj4string(stack.targets) ) ) ) {
      stop(sprintf("%s: projections of maps stack and training/validation shapefiles do not match.", zone))
    }

    poly.validation <- poly.training.validation[poly.training.validation@data$TV == "VALIDATION", ]
    
    ##-----------
    # poly.validation <- poly.validation[1:40, ]
    ##-----------
    
    sizeBlocks <- ceiling( nrow(poly.validation@data) / min(detectCores(), nrow(poly.validation@data)) )  ## to avoid problems when number of rows of the matrix is smaller than the nr of cores 
    nblocks <- ceiling(nrow(poly.validation@data)/sizeBlocks)   ## to avoid problems when sizeBlocks is small (nblocks will allow to cover just a little more, or exactly, the size of data)
    nr.clusters <- nblocks
    cl <- makeCluster(nr.clusters)
    registerDoParallel(cl)
    targets.extract <- foreach (blck = 1:nblocks, .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=T for more info when debugging
      ## WHEN NOT USING FOREACH, TESTING PHASE
    # nblocks <- 19
    # for (blck in 17:nblocks) {
    #   sizeBlocks <- 20
      ## WHEN NOT USING FOREACH, TESTING PHASE
      
      indPart <- seq(from=(blck-1)*sizeBlocks+1, to=min(blck*sizeBlocks, nrow(poly.validation@data)), by=1)
  
      sampling.shp <- poly.validation[indPart, ]   ## for actual study we use polygons
      extract(stack.targets, sampling.shp, fun = mean, na.rm = F, weights = T, normalizeWeights = T,
                            cellnumbers = F, small = F, df = T, factors = F, sp = F)  # extract is a function that actually accesses the values of the rasters, MEMORY ISSUES?
    }
    stopCluster(cl)
    
    ## check if nr of rows is the same in input and output of the foreach (to avoid issues with block indexes)
    if ( nrow(poly.validation@data) != nrow(targets.extract) ) {
      stop(sprintf("Zone %s: nr. of rows is different for poly.validation and targets.extract", zone))
    }
    
    targets.extract$ID <- poly.validation@data$POLY250ID
    
    full.table.targets[[z]] <- targets.extract
    
    ## clock UTM zone time
    temp.toc <- proc.time()-temp.tic[3]
    print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))
  
  }  ## end if on params4$check.plot.vs.pix
  

#### EXTRACT ELEV_P95 TEMPORAL TRAJECTORIES ON VAL --------------------------------------------
  
  if (params4$check.temporal.val.traj) { 
    
    ## read 1984-2012 maps for elev_p95
    for (yr in params4$mapped.years) {
      cmd <- sprintf('map%s <- raster(file.path("%s", "UTM_results", "UTM_%s_elev_p95_RF_1984_2012", "UTM_%s_elev_p95_RF_%s.dat", fsep = .Platform$file.sep))', yr, Landsat_dir, zone.nr, zone.nr, yr)  
      eval(parse(text=cmd))
    }
    
    vars.years <- c(params4$mapped.years)  
    str.years.vect <- sprintf("map%s", params4$mapped.years)
    str.years <- paste(str.years.vect, collapse = ', ')
    
    cmd <- sprintf("stack.years <- stack(%s)", str.years)  
    eval(parse(text=cmd))
    names(stack.years) <- vars.years
    
    ## check if coordinate systems are the same between raster stack and shapefiles
    if ( !all( sapply( list(proj4string(poly.training.validation)), FUN = identical, proj4string(stack.years) ) ) ) {
      stop(sprintf("%s: projections of yearly maps stack and training/validation shapefiles do not match.", zone))
    }
    
    poly.validation <- poly.training.validation[poly.training.validation@data$TV == "VALIDATION", ]
    cpt.poly.validation <- cpt.poly.training.validation[cpt.poly.training.validation@data$TV == "VALIDATION", ]
    
    ##-----------
    # poly.validation <- poly.validation[1:40, ]
    # cpt.poly.validation <- cpt.poly.validation[1:40, ]
    ##-----------
    
    sizeBlocks <- ceiling( nrow(poly.validation@data) / min(detectCores(), nrow(poly.validation@data)) )  ## to avoid problems when number of rows of the matrix is smaller than the nr of cores 
    nblocks <- ceiling(nrow(poly.validation@data)/sizeBlocks)   ## to avoid problems when sizeBlocks is small (nblocks will allow to cover just a little more, or exactly, the size of data)
    nr.clusters <- nblocks
    cl <- makeCluster(nr.clusters)
    registerDoParallel(cl)
    years.extract <- foreach (blck = 1:nblocks, .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=T for more info when debugging
      ## WHEN NOT USING FOREACH, TESTING PHASE
      # nblocks <- 19
      # for (blck in 17:nblocks) {
      #   sizeBlocks <- 20
      ## WHEN NOT USING FOREACH, TESTING PHASE
      
      indPart <- seq(from=(blck-1)*sizeBlocks+1, to=min(blck*sizeBlocks, nrow(poly.validation@data)), by=1)
      
      ## selection of the sampling unit
      if (params4$temp.sampling.shp == 'poly') {   
        sampling.shp <- poly.validation[indPart, ]   ## for actual study we use polygons
      } else if (params4$temp.sampling.shp == 'cpt') { 
        sampling.shp <- cpt.poly.validation[indPart, ]    ## to speed up tests read only centerpoint
      }
      
      extract(stack.years, sampling.shp, fun = mean, na.rm = F, weights = T, normalizeWeights = T,
              cellnumbers = F, small = F, df = T, factors = F, sp = F)  # extract is a function that actually accesses the values of the rasters, MEMORY ISSUES?
    }
    stopCluster(cl)
    
    ## check if nr of rows is the same in input and output of the foreach (to avoid issues with block indexes)
    if ( nrow(poly.validation@data) != nrow(years.extract) ) {
      stop(sprintf("Zone %s: nr. of rows is different for poly.validation and years.extract", zone))
    }
    
    years.extract$ID <- poly.validation@data$POLY250ID  ## identical to cpt.poly.validation@data$POLY250ID
    
    full.table.years[[z]] <- years.extract
    
  }  ## end if on params4$check.temporal.val.traj | params4$check.temporal.val.boxpl
  
}   ## end for on paramsGL$zones

## subdirectory to save model assessment plots
MapCheck.CAN.subdir <- file.path(base_figures_dir, "MapCheck_CAN_level", fsep = .Platform$file.sep)
if (! file.exists(MapCheck.CAN.subdir)){dir.create(MapCheck.CAN.subdir, showWarnings = F, recursive = T)}

#### 1) PIXEL VS PLOT SCATTERPLOT ---------------------------------------------------

if (params4$check.plot.vs.pix) {
  
  ## subdirectory to save plot vs pixel plots
  plot.vs.pix.CAN.subdir <- file.path(MapCheck.CAN.subdir, "Plot_VS_Pixel", fsep = .Platform$file.sep)
  
  if (! file.exists(plot.vs.pix.CAN.subdir)){dir.create(plot.vs.pix.CAN.subdir, showWarnings = F, recursive = T)}
  
  target.predictions <- as.data.table(rbindlist(full.table.targets))
  colnames(target.predictions)[1] <- "FCID"
  
  target.predictions.contin <- maps_int_2_contin(target.predictions[, 2:ncol(target.predictions), with=FALSE], params4$multipl) 
  target.predictions.contin[, FCID:=target.predictions$FCID]
  
  dt.to.compare <- merge(target.predictions.contin, plot.based.predictions@data, by="FCID")
  
  idx.targ <- 1
  for (targ in params3$targ.names.sh) {
    
    cmd <- sprintf('x.axis.var <- dt.to.compare$P_%s', params3$targ.names.sh[idx.targ])  ## values predicted on plots ("P_") on the x axis
    eval(parse(text=cmd))
    cmd <- sprintf('y.axis.var <- dt.to.compare$%s', params3$targ.names.sh[idx.targ])  # values predicted on pixels on y axis
    eval(parse(text=cmd))
  
    ## CAN-level scatterplot of plot vs pixel based predictions
    scatt.range <- c(min(x.axis.var, y.axis.var), quantile(c(x.axis.var, y.axis.var), 1, names=F))
    fig.name.str <- file.path(plot.vs.pix.CAN.subdir, sprintf("PlotPixel_scatter_%s.pdf", targ), sep='')
    title.str <- sprintf("%s [%s]: R^2=%.3f", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ], cor(x.axis.var, y.axis.var)**2)
    pdf(fig.name.str)
      theme_set(theme_gray(base_size = 18))
      scatt.plot <- plot_colorByDensity(x.axis.var, y.axis.var, xlim=scatt.range, ylim=scatt.range, xlab="plot prediction", ylab="pixel prediction", main=title.str) +
        geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line
      print(scatt.plot)
    dev.off()
  
    idx.targ <- idx.targ + 1
  
  }
  
  write.csv(arrange(dt.to.compare, FCID), file = file.path(plot.vs.pix.CAN.subdir, "plot_vs_pixel_predictions.csv", sep = ''))

}  ## end if on params4$check.plot.vs.pix

if (params4$check.temporal.val.traj) {
  
  ## load explanatory and response variables extracted for training plots
  X.trn.val.raw <- read.csv(file.path(base_wkg_dir, "X_trn_val.csv", fsep = .Platform$file.sep))
  Y.trn.val.raw <- read.csv(file.path(base_wkg_dir, "Y_trn_val.csv", fsep = .Platform$file.sep))
  
  yearly.predictions <- as.data.frame(rbindlist(full.table.years))
  colnames(yearly.predictions)[1] <- "FCID"
  yearly.predictions$FCID <- as.character(yearly.predictions$FCID) 
  X.trn.val.raw$FCID <- as.character(X.trn.val.raw$FCID) 
  
  ## merge the three dataframes by FCID to have corresponding values of pixel and plot based predictions
  X.trn.val.raw.w.ECO <- inner_join(x = X.trn.val.raw, y = Y.trn.val.raw[, c("FCID", "ECOZONE")], by="FCID")
  X.yearly.predictions <- inner_join(x = X.trn.val.raw.w.ECO, y = yearly.predictions, by="FCID")
  df.with.year <- inner_join(x = plot.based.predictions@data[, c("FCID", "O_el_p95")], y = X.yearly.predictions, by="FCID")
  
#### 2) TEMPORAL TRAJECTORIES ON VAL SET -----------------------------------------------
  
  if (params4$check.temporal.val.traj) {
  
    ## summarize all columns (where possible, otherwise NA) with mean
    res <- df.with.year %>% 
           group_by(YrsSince_GrCh, Ch_attr) %>% 
           summarise_each(funs(mean))
    ## add counts for each combination of YrsSince_GrCh and Ch_attr
    res2 <- df.with.year %>% 
            group_by(YrsSince_GrCh, Ch_attr) %>% 
            summarise(n = n())
    mean.by.yrCh.Chattr <- data.frame(count=res2$n, res) %>% ## put them together while keeping only counts higher than 8
                           filter(count >= 8)
    
    temporal.traj <- mean.by.yrCh.Chattr[, as.character(paste("X", params4$mapped.years,  sep=''))]
    colnames(temporal.traj) <- gsub("X", "", colnames(temporal.traj) )
    
    year.of.change <- paramsGL$TARGET_YEAR-mean.by.yrCh.Chattr$YrsSince_GrCh
    temporal.traj.melted <- melt(temporal.traj)
    temporal.traj.melted$rowid <- sprintf("%s in %s", params3$Ch_attr.classes[mean.by.yrCh.Chattr$Ch_attr+1], year.of.change)
    # temporal.traj.melted <- melt(temporal.traj[1:10, ])
    # temporal.traj.melted$rowid <- 1:10
    
    y.to.highlight <- rep(NA, nrow(mean.by.yrCh.Chattr))
    i <- 1:nrow(mean.by.yrCh.Chattr)
    j <- match( year.of.change, gsub("X", "", colnames(mean.by.yrCh.Chattr)) )
    y.to.highlight <- mean.by.yrCh.Chattr[cbind(i, j)]
  
    idx.targ <- which(params3$targ.names.sh == params4$targ.to.plot)
    
    ## temporal trajectories by Yr_Gr_Change and Change_attr
    fig.name.str <- file.path(MapCheck.CAN.subdir, sprintf("TempTraj_elev_p95.pdf"), sep='')
    # title.str <- sprintf("%s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ])
    temporal.traj.melted$variable <- as.numeric(as.character(temporal.traj.melted$variable))
    scale.ten.classes <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
    pdf(fig.name.str)
    theme_set(theme_gray(base_size = 14))
    time.series <- 
      ggplot(data=NULL) + 
      geom_line(data=temporal.traj.melted, aes( variable, value, size = 4, group=factor(rowid), color=factor(rowid) ) ) +
      scale_color_manual(values=scale.ten.classes) +
      scale_x_continuous(breaks=seq(1984, 2012, by=2), lim=c(1984,2012)) +
      xlab("year") +
      ylab(sprintf("%s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ])) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_point( data = data.frame(x=year.of.change, y=y.to.highlight), aes(x=x, y=y, size = 5, color=factor(unique(temporal.traj.melted$rowid))) ) +
      scale_color_manual(values=scale.ten.classes) 
    print(time.series)
    dev.off()
  
  }   ## end if on params4$check.temporal.val.traj
  
}  ## end if on params4$check.temporal.val.traj | params4$check.temporal.val.boxpl

#### 3) MAP COHERENCE -------------------------------------------------------

if (params4$check.map.coherence) {
  
  # removeTmpFiles(h = 0)   ## remove all temporary files found in the temp directory ( dirname(rasterTmpFile()) ) 
  
  # scaled.el_p95 <- scale(stack.targets$el_p95)
  # scaled.p_1r_2m <- scale(stack.targets$p_1r_2m)
  # incoherence <- scaled.el_p95-scaled.p_1r_2m
  # lims.el_p95 <- quantile(stack.targets$el_p95, probs=params4$incoher.quant)
  # lims.p_1r_2m <- quantile(stack.targets$p_1r_2m, probs=params4$incoher.quant)
  # stats.scaled <- c(cellStats(scaled.el_p95, 'mean'), cellStats(scaled.el_p95, 'sd'), cellStats(scaled.p_1r_2m, 'mean'), cellStats(scaled.p_1r_2m, 'sd'))
  
  ht.low.cover.high <- stack.targets$el_p95 < params4$lims.el_p95[1] & stack.targets$p_1r_2m > params4$lims.p_1r_2m[2]
  ht.high.cover.low <- stack.targets$el_p95 > params4$lims.el_p95[2] & stack.targets$p_1r_2m < params4$lims.p_1r_2m[1]
  
  # writeRaster(incoherence, file.path(MapCheck.CAN.subdir, sprintf("Standard_elev_p95_minus_p_1r_2m"), sep=''), format='GTiff', overwrite=T)
  writeRaster(ht.low.cover.high, file.path(MapCheck.CAN.subdir, "Maps", sprintf("height_below_%s_cover_above_%s", params4$lims.el_p95[1], params4$lims.p_1r_2m[2]), sep=''), format='GTiff', datatype="LOG1S", overwrite=T)
  writeRaster(ht.high.cover.low, file.path(MapCheck.CAN.subdir, "Maps", sprintf("height_above_%s_cover_below_%s", params4$lims.el_p95[2], params4$lims.p_1r_2m[1]), sep=''), format='GTiff', datatype="LOG1S", overwrite=T)
  
}  ## end if on params4$check.map.coherence


## delete temp folder with temp rasters
unlink(raster_tmp_dir, recursive = T, force = T)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog4, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
