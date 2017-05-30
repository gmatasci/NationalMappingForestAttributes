## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog6_temporal_analysis.R                             
## Objective: Plot time-series of mapped attributes with boxplots of an ensemble of pixels sampled over large areas or with lines for individual pixel

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# If params6$check.temporal.val <- T: 
# - use foreach on blocks
# - check no redefinition of paramsGL$zones and subsetting of poly.validation: look for ##---
# - check parameters

# - check influence of rasterOptions(maxmemory = 1e+09), rasterOptions(chunksize = 1e+08)


## SOLVED:
# -V change back and remove _v2 -- now we have separate folders for each version of the csvs
# -V save separate correlation trends plots -- plots are now produced
# -V Check for 0,0,0,0 in bands values: remove them! -- removed in prog4c with line  yearly.predictors <- yearly.predictors[b1!=0 & b2!=0 & b3!=0 & b4!=0 & b5!=0 & b6!=0]
# -V fix pixID issue band5 wrt rest of plots -- order of sampled pixels is prreserved by using merge to keep matching pixel IDs instead of cbinding data.tables which can lead to mixing different orders
# -V implement and check script for new boxplots v3 -- v3 boxplots are more stable bc of the further filtering of the sampled pixels


#### READS/WRITES ------------------------------------------------------------


#### INIT --------------------------------------------------------------------

print('Prog6: temporal analysis')

rm(list=ls()) ## clear all variables

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params2.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params3.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params6a.Rdata"
load(param_file)


source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")


#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params6 <- list()

params6$check.temporal.val <- F  ## wheter to run the analysis over time on val set
params6$temporal.val.data.source <- "CSV"  ## whether to read val set values from CSVs produced by sampling with IDL scripts (option "CSV") or from the raw maps in R using the shapefiles (option "MAP")
params6$temp.val.sampling.shp <- "cpt"    ## sampling unit (to be set as "poly") to read back map values on validation plots, to be set as "cpt" to speed things up in development 
paramsGL$zones <- "UTM13S"   ## redefine UTM zones to consider for params6$temporal.val.data.source <- "MAP" case

params6$check.temporal.full <- T  ## wheter to run the analysis over time on full maps
# params6$temporal.data.types <- c("raw", "fitted")  ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends
params6$temporal.data.types <- c("fitted")

# params6$change.types <-  c("always_treed", "fire", "harvesting")
params6$change.types <- c("fire", "harvesting")
# params6$change.types <-  c("harvesting")

params6$methods <- c("YAI")   ## for which methods to run the analysis, can accept both c("RF", "YAI") 

params6$mapped.years <- seq(from=1984, to=2012)   ## only used for val boxplots when reading back the images
# params6$mapped.years <- seq(from=1984, to=1988)
# params6$mapped.years <- seq(from=2008, to=2012)

params6$forest.classes <- list("AllForest"=c(81, 210, 220, 230), "WetlandTreed"=81, "ConiferousForest"=210, "DeciduousForest"=220, "MixedForest"=230)
# params6$forest.classes <- list("AllForest"=c(81, 210, 220, 230), "ConiferousForest"=210)

# params6$start.change.plots <- -50
params6$start.change.plots <- -3

params6$plot.count.hist <- T   ## whether to produce histogram of counts in each annual bin in the boxplots time-series

params6$plot.smoothed.ts <- F    ## whether to produce Loess smoothed plots of single pixel time-series

# params6$mapped.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West",
#                                  "Hudson Plains", "Montane Cordillera", "Pacific Maritime", "Taiga Cordillera",
#                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")
# params6$ecozone.codes <- c(2, 3, 4, 5, 6, 7, 9, 11, 15, 16, 17, 18)
params6$mapped.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", "Montane Cordillera", "Pacific Maritime")
params6$ecozone.codes <- c(2, 3, 4, 5, 6, 9, 11)
# params6$mapped.ecozones <- list("Pacific Maritime")
# params6$ecozone.codes <- c(11)

params6$targ.names <- c("elev_p95", "percentage_first_returns_above_2m")
params6$targ.units <- c("m", "%")
params6$targ.names.plots <- c("elev_p95", "cover_2m")
params6$ylim.targ$absolute$elev_p95$NONBOREAL <- c(0, 50)
params6$ylim.targ$absolute$elev_p95$BOREAL <- c(0, 28)
params6$ylim.targ$absolute$cover_2m <- c(0, 100)
params6$ylim.targ$difference$elev_p95 <- c(-11, 11)
params6$ylim.targ$difference$cover_2m <- c(-40, 40)
params6$ylim.band5 <- c(200, 3000)

params6$nr.pix.time.series <- 12

param_file_prog6 = file.path(base_wkg_dir, 'AllUTMzones_params6.Rdata', fsep = .Platform$file.sep)
save(params6, file = param_file_prog6)

#### LOAD PACKAGES ----------------------------------------------------------

## install data.table version with fwrite
# install.packages("data.table", repos = "https://Rdatatable.github.io/data.table", type = "source")

list.of.packages <- c("ggplot2",
                      "gridExtra",
                      "grid",
                      "reshape2",
                      "rgdal",
                      "raster",
                      "GSIF",
                      "sp",
                      "dplyr",    ## to be loaded before foreach to avoid "assertion failed" errors
                      "plyr", 
                      "landsat",
                      "lattice",
                      "lubridate", 
                      "doParallel",
                      "foreach",
                      "data.table",
                      "methods"
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

#### 1) BOXPLOTS ON VAL SET ---------------------------------------------------

## Run only if we check for temporal boxplots VS 2010 observed values
if (params6$check.temporal.val) {  

  print('Validation set boxplots:')
  
  ## Update results directory to point to the FINAL results one
  base_results_dir <- file.path(base_dir, "FINAL_RUN", "results", fsep=.Platform$file.sep)
  
  ValBoxPlots.subdir <- file.path(Temporal.subdir, "ValBoxPlots", fsep = .Platform$file.sep)
  if (! file.exists(ValBoxPlots.subdir)){dir.create(ValBoxPlots.subdir, showWarnings = F, recursive = T)}
  
  val.time.series.folder.name <- "E:/NTEMS/UTM_results/samples_plots/predicted_values"
  
  ## Load predictions on Val set 
  plot.based.predictions <- readOGR(dsn = base_results_dir, layer = 'val_predictions_RF')
  colnames(plot.based.predictions@data)[ncol(plot.based.predictions@data)] <- "FCID"
  plot.based.predictions@data$FCID <- as.character(plot.based.predictions@data$FCID) 

  ## Load explanatory and response variables extracted for training plots
  X.trn.val.raw <- fread(file.path(base_wkg_dir, "X_trn_val.csv", fsep = .Platform$file.sep))
  Y.trn.val.raw <- fread(file.path(base_wkg_dir, "Y_trn_val.csv", fsep = .Platform$file.sep))
  
  ## Merge by FCID
  X.trn.val.raw$FCID <- as.character(X.trn.val.raw$FCID) 
  X.trn.val.raw.w.ECO <- merge(X.trn.val.raw, Y.trn.val.raw[, c("FCID", "ECOZONE"), with=FALSE], by="FCID")
  
  ## Load infos about plots
  plots.info <- fread("D:/Research/ANALYSES/NationalMappingForestAttributes/Rcode_supercomputer_Txomin/Get_plots_coords/plot_UTM_info.csv", header = T)
  
  if (params6$temporal.val.data.source == "MAP") {
  
    full.table.targets <- vector("list", length(paramsGL$zones)) 
    full.table.years <- vector("list", length(paramsGL$zones)) 
    for (z in 1:length(paramsGL$zones)) {
      
      temp.tic <- proc.time() ## start clocking time for each UTM zone
      
      zone <- paramsGL$zones[z]
      zone.nr <- substr(zone, 4, nchar(zone))
      print(paste('Processing zone', zone))   ## converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
      wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)
    
      ## Read training and validation shapefiles: polygons and centerpoints
      poly.training.validation <- readOGR(dsn = wkg_dir, layer = paste(zone,"_poly_training_validation",sep = ''))
      cpt.poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_pt_centerpt_training_validation",sep = ''))

#### EXTRACT ELEV_P95 TEMPORAL TRAJECTORIES ON VAL --------------------------------------------
      
      ## Read 1984-2012 maps for elev_p95
      for (yr in params6$mapped.years) {
        cmd <- sprintf('map%s <- raster(file.path("%s", "UTM_results", "RF_corrected", "UTM_%s_elev_p95_RF_%s.dat", fsep = .Platform$file.sep))', yr, Landsat_dir, zone.nr, yr)  
        eval(parse(text=cmd))
      }
      
      vars.years <- c(params6$mapped.years)  
      str.years.vect <- sprintf("map%s", params6$mapped.years)
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
        
        ## Selection of the sampling unit
        if (params6$temp.val.sampling.shp == 'poly') {   
          sampling.shp <- poly.validation[indPart, ]   ## for actual study we use polygons
        } else if (params6$temp.val.sampling.shp == 'cpt') { 
          sampling.shp <- cpt.poly.validation[indPart, ]    ## to speed up tests read only centerpoint
        }
        
        extract(stack.years, sampling.shp, fun = mean, na.rm = F, weights = T, normalizeWeights = T,
                cellnumbers = F, small = F, df = T, factors = F, sp = F)  # extract is a function that actually accesses the values of the rasters, MEMORY ISSUES?
      }
      stopCluster(cl)
      
      ## Check if nr of rows is the same in input and output of the foreach (to avoid issues with block indexes)
      if ( nrow(poly.validation@data) != nrow(years.extract) ) {
        stop(sprintf("Zone %s: nr. of rows is different for poly.validation and years.extract", zone))
      }
      
      years.extract$ID <- poly.validation@data$POLY250ID  ## identical to cpt.poly.validation@data$POLY250ID
      
      full.table.years[[z]] <- years.extract
        
    }   ## end for on paramsGL$zones
  
    yearly.predictions.R <- as.data.table(rbindlist(full.table.years))
    
    colnames(yearly.predictions.R)[1] <- "FCID"
    yearly.predictions.R$FCID <- as.character(yearly.predictions.R$FCID) 
    
    setnames(yearly.predictions.R, old=sprintf("X%s", params6$mapped.years), new=as.character(params6$mapped.years))
    yearly.predictions.R <- data.table(FCID=yearly.predictions.R[, FCID], yearly.predictions.R[, as.character(params6$mapped.years), with=FALSE]/1000)  ## convert back to proper values
    
  }  ## end if params6$temporal.val.data.source == "MAP"
  
  
  
#### PLOT BOXPLOTS ---------------------------------------------------------------------
  
  for (targ in params6$targ.names) {
    
    print(sprintf("  Targ: %s", targ))
    
    idx.targ <- which(targ == params6$targ.names)   ## numeric index useful to select elements of string vectors
    
    ## set of target-specific parameters
    if (targ == "elev_p95") {
      observed.val.col <- "O_el_p95"  ## name of the column with the observed values in plot.based.predictions@data
      ylim.targ.abs <- params6$ylim.targ$absolute$elev_p95
      ylim.targ.diff <- params6$ylim.targ$difference$elev_p95
    } else if (targ == "percentage_first_returns_above_2m") {
      observed.val.col <- "O_p_1r_2m"
      ylim.targ.abs <- params6$ylim.targ$absolute$cover_2m
      ylim.targ.diff <- params6$ylim.targ$difference$cover_2m
    }
    
    for (method in params6$methods) { 
      
      print(sprintf("    Method: %s", method))
      
      print("      Ch_attr: No change")
      
      if (params6$temporal.val.data.source == "CSV") {
        file.name <- file.path(val.time.series.folder.name, sprintf("time_series_val_set_%s_%s.csv", method, targ), fsep = .Platform$file.sep)
        yearly.predictions.CSV <- fread(file.name, header = T)
        yearly.predictions.CSV[, FCID:=plots.info$FCID[plots.info$TV=="VALIDATION"]]
        yearly.predictions <- yearly.predictions.CSV
      } else if (params6$temporal.val.data.source == "MAP") {
        if (method != "RF" | targ != "elev_p95") {next}
        yearly.predictions <- yearly.predictions.R
      }
      
      ## Check match between plot and pixel (no weigthed average) based predictions
      if (method == "RF" & targ == "elev_p95") {

        ## On all Val set
        pix.pred <- yearly.predictions[, c("FCID", "2010"), with=FALSE]
        plot.pred <- data.table(PlotPred=plot.based.predictions@data$P_el_p95, FCID=plot.based.predictions$FCID)
        merged.pred <- merge(pix.pred, plot.pred, by="FCID")
        setnames(merged.pred, old="2010", new="PixelPred")
        R2coeff <- cor(merged.pred$PlotPred, merged.pred$PixelPred)**2
        fig.name.str <- file.path(base_figures_dir, "MapCheck_CAN_level", "Plot_VS_Pixel_CAN", sprintf("Plot_VS_Pixel_CAN_no_weigh_avg_%s_%s.pdf", params6$targ.names.plots[idx.targ], method), sep='')
        pdf(fig.name.str)
           p <- plot_colorByDensity(merged.pred$PlotPred, merged.pred$PixelPred, xlab="plot prediction", ylab="pixel prediction", main=sprintf("%s [%s]: R^2=%.3f on %s validation points", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ], R2coeff, nrow(merged.pred))) +
                geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line                        
        print(p)
        dev.off()
        
        ## On UTM zones 12S and 13S only (to compare with reading back preliminary maps)
        subset.UTM12S13S <- X.trn.val.raw.w.ECO %>%
          filter(TV=="VALIDATION" & (UTMzone=="13S" | UTMzone=="12S") )
        subset.UTM12S13S <- as.data.table(subset.UTM12S13S)
        merged.pred.subs <- merge(subset.UTM12S13S, merged.pred, by="FCID")
        R2coeff <- cor(merged.pred.subs$PlotPred, merged.pred.subs$PixelPred)**2
        fig.name.str <- file.path(base_figures_dir, "MapCheck_CAN_level", "Plot_VS_Pixel_UTM12S13S", sprintf("Plot_VS_Pixel_UTM12S13S_no_weigh_avg_%s_%s.pdf", params6$targ.names.plots[idx.targ], method), sep='')
        pdf(fig.name.str)
           p <- plot_colorByDensity(merged.pred.subs$PlotPred, merged.pred.subs$PixelPred, xlab="plot prediction", ylab="pixel prediction", main=sprintf("%s [%s]: R^2=%.3f on %s validation points", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ], R2coeff, nrow(merged.pred.subs))) +
                geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line
        print(p)
        dev.off()
        
      }
      
      ## merge the three dataframes by FCID to have corresponding values of pixel and plot based predictions
      X.yearly.predictions <- merge(X.trn.val.raw.w.ECO, yearly.predictions, by="FCID")
      dt.with.year <- merge(as.data.table(plot.based.predictions@data[, c("FCID", observed.val.col)]), X.yearly.predictions, by="FCID")
      
      ## Predicted boxplots over time on validation set + observed boxplot in 2010
      for (ez in 1:length(params6$mapped.ecozones)) {
        
        ecozone.name <- params6$mapped.ecozones[ez]
        
        ## Skip if the ecozone in not among those found in the plots 
        if ( !ecozone.name %in% unique(dt.with.year$ECOZONE) ) { next } 
        
        ## Select only no-change plots of the ecozone of interest
        dt.final <- dt.with.year %>%
                    filter(ECOZONE==ecozone.name, YrsSince_GrCh==50, Ch_attr==0)
        dt.final <- as.data.table(dt.final)
        dt.final <- dt.final[, c(params6$mapped.years, observed.val.col), with=FALSE]
        
        ## Rename columns to reflect predicted (P_) and observed (O_) values
        colnames(dt.final) <- c( sprintf("P_%s", colnames(dt.final)[1:length(colnames(dt.final))-1]), "O_2010")
        
        ## Computing residuals wrt to 2010 observed value
        dt.final.point.diff <- as.data.table(apply(dt.final, 2, FUN=function (x) x-dt.final$O_2010 ))
        
        yearly.predictions <- melt(dt.final)
        yearly.predictions.point.diff <- melt(dt.final.point.diff)
        
        ## Boxplot of predicted values compared to observed values in 2010 
        fig.name.str <- file.path(ValBoxPlots.subdir, sprintf("NoChange_ValSet_%s_%s_%s_1984_2012.pdf", gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
        title.str <- sprintf("Val. set:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, nrow(dt.final))
        pdf(fig.name.str)
        theme_set(theme_gray(base_size = 11))
        p <- ggplot(yearly.predictions) +
          geom_boxplot(aes(x=variable, y=value), outlier.shape=NA, notch=F) +
          ggtitle(title.str) +
          coord_cartesian(ylim=ylim.targ.abs) +   ## don't use ylim() as it will change shape bc it remove points beyond limit !!!
          xlab("year") +
          ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
        print(p)
        dev.off()
        
        # Difference with observed 2010 value (residuals wrt 2010)
        fig.name.str <- file.path(ValBoxPlots.subdir, sprintf("NoChange_ValSet_Residuals_%s_%s_%s_1984_2012.pdf", gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
        title.str <- sprintf("Per-sample difference (residual) wrt 2010 \n Val. set:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, nrow(dt.final))
        pdf(fig.name.str)
        theme_set(theme_gray(base_size = 11))
        p <- ggplot(yearly.predictions.point.diff) +
          geom_boxplot(aes(x=variable, y=value), outlier.shape=NA, notch=F) +
          ggtitle(title.str) +
          coord_cartesian(ylim=ylim.targ.diff) +   ## don't use ylim() as it will change shape bc it remove points beyond limit !!!
          xlab("year") +
          ylab(sprintf("residual %s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
        print(p)
        dev.off()
          
      }   ## end for on params6$ecozone.code
      
    }  ## end for on params6$methods
     
  }   ## end for on params6$targ.names
  
}  ## end if params6$check.temporal.val



#### 2) PLOTS ON FULL MAP ---------------------------------------------------
  
if (params6$check.temporal.full) {
  
  print('Full map plots:')
  
  for (temporal.data.type in params6$temporal.data.types) {
    
    print(sprintf("Temporal data: %s", temporal.data.type))
    
    for (change.type in params6$change.types) {
      
      print(sprintf("    Ch_attr: %s", change.type))
      
      samples.dir <- file.path(params6a$temporal.base.dir, change.type, 'samples_u', fsep = .Platform$file.sep)
      landcover.dir <- file.path(params6a$temporal.base.dir, change.type, 'landcover', fsep = .Platform$file.sep)
      
      if (temporal.data.type == "fitted") {
        predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics_fitted', fsep = .Platform$file.sep)
        predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, 'predicted_values_fitted_NOTEMPORAL', fsep = .Platform$file.sep)
        Temporal.subdir <- file.path(sprintf('%s_FITTED_NOTEMPORAL', base_figures_dir), "TemporalAnalyses", fsep = .Platform$file.sep)
      } else if (temporal.data.type == "raw") {
        predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics', fsep = .Platform$file.sep) 
        predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, 'predicted_values', fsep = .Platform$file.sep) 
        Temporal.subdir <- file.path(base_figures_dir, "TemporalAnalyses", fsep = .Platform$file.sep)
      }
      FullBoxPlots.subdir <- file.path(Temporal.subdir, "FullBoxPlots", change.type, fsep = .Platform$file.sep)
      FullPixelTimeSeries.subdir <- file.path(Temporal.subdir, "FullPixelTimeSeries", change.type, fsep = .Platform$file.sep)
      if (! file.exists(FullBoxPlots.subdir)){dir.create(FullBoxPlots.subdir, showWarnings = F, recursive = T)}
      if (! file.exists(FullPixelTimeSeries.subdir)){dir.create(FullPixelTimeSeries.subdir, showWarnings = F, recursive = T)}
      
      for (ez in 1:length(params6$mapped.ecozones)) {
        
        ecozone.name <- params6$mapped.ecozones[ez]
        ecozone.code <- params6$ecozone.codes[ez]
        
        print(sprintf("        Ecozone: %s", ecozone.name))
        
        for (targ in params6$targ.names) {
            
          print(sprintf("            Targ: %s", targ))
          
          idx.targ <- which(targ == params6$targ.names)   ## numeric index useful to select elements of string vectors

          ## set of target and ecozone-specific parameters
          if (targ == "elev_p95") {
            if (ecozone.name %in% params3$sampled.ecozones.NONBOREAL) {
              ylim.targ.abs <- params6$ylim.targ$absolute$elev_p95$NONBOREAL
            } else {
              ylim.targ.abs <- params6$ylim.targ$absolute$elev_p95$BOREAL
            }
          } else if (targ == "percentage_first_returns_above_2m") {
            ylim.targ.abs <- params6$ylim.targ$absolute$cover_2m
          }
          
          for (method in params6$methods) { 
            
            print(sprintf("                Method: %s", method))
    
            file.name.melted.predictions <- file.path(predicted.values.dir, sprintf("%s_%s_%s_%s.csv", change.type, method, targ, gsub(' ', '', ecozone.name)), fsep = .Platform$file.sep)
            yearly.predictions.full <- fread(file.name.melted.predictions, header = T)
            yearly.predictions.full <- yearly.predictions.full[, -"ecozone", with=FALSE]
            setkey(yearly.predictions.full, pixID)
            
#### MERGE AND PLOT BAND 5 DATA ---------------------------------------------------
            
            ## Build band 5 data time-series for comparison with time-series of predicted data and get year of change info and land-cover class for each pixel
            if (targ == params6$targ.names[1] & method == params6$methods[1]) {  ## to do this plot only once per ecozone
              
              file.name <- file.path(samples.dir, sprintf("%s_%s.csv", ecozone.code, ecozone.name), fsep = .Platform$file.sep)
              pix.IDs <- fread(file.name, header = F)
              
              file.name <- file.path(landcover.dir, sprintf("%s_%s_vlce.csv", ecozone.code, ecozone.name), fsep = .Platform$file.sep)
              landcover <- fread(file.name, header = F)
              landcover[, pixID:=sprintf('%s_%s', landcover[,1][[1]], landcover[,2][[1]])]
              landcover[,1:2] <- NULL
              colnames(landcover)[1:ncol(landcover)-1] <- params6$mapped.years
              setkey(landcover, pixID)
              
              ## Loop over the years to merge band 5 values of each year 
              band.year.dt <- data.table(pixID=yearly.predictions.full$pixID, order=seq(from=1, to=nrow(yearly.predictions.full)))  ## initialize the containing dt with pixIDs from table with predicted values read from predicted.values.dir (can be smaller than sample info table read from samples.dir) 
              setkey(band.year.dt, pixID)
              for (yr in 1:length(params6$mapped.years)) {
                
                year <- params6$mapped.years[yr]
                
                ## Load dataset with bands info from CSV (predictors dataset)
                if (temporal.data.type == 'fitted') {
                  file.name <- file.path(predictors.dir, sprintf("%s_%s_%s_fitted.csv", ecozone.code, ecozone.name, year), fsep = .Platform$file.sep)
                } else {
                  file.name <- file.path(predictors.dir, sprintf("%s_%s_%s.csv", ecozone.code, ecozone.name, year), fsep = .Platform$file.sep)
                }
                
                yearly.predictors <- fread(file.name, header = T)  ## same size and order of samples as pix.IDs (sample info)
                yearly.predictors[, pixID:=sprintf('%s_%d', pix.IDs[[1]], pix.IDs[[2]])]  ## so we can attach unique pixel ID with UTM zone info and index from pix.IDs
                yearly.predictors[, index:=NULL]
                setkey(yearly.predictors, pixID)  ## set key for merges as pixID
                
                ## Merge band 5 values based on pixID
                band.year.dt <- merge(band.year.dt, yearly.predictors[, .(pixID, b5)], all=FALSE)
                colnames(band.year.dt)[ncol(band.year.dt)] <- as.character(year)
                
                ## Attach year of change info to band.year.dt, as this is the dt that stays unchanged irrespective of method and targ 
                if (yr == length(params6$mapped.years)) {   ## do it on last year, as only then we know with certainty when a pixel has changed
                  if (year != 2012) {
                    warning("Not using last year to determine YrOfChange")
                  }
                  band.year.dt <- merge(band.year.dt, yearly.predictors[, .(pixID, YrsSince_GrCh)], all=FALSE)
                  band.year.dt[, YrOfChange:= year - YrsSince_GrCh]
                }
                
              }   ## end for on params6$mapped.years
              
              setkey(band.year.dt, order)  ## reset order
              band.year.dt[, "order":=NULL]   ## ...and then remove it
  
              ## Get pixel indices to plot in time-series
              set.seed(paramsGL$global.seed)
              pix.to.plot.idx <- sample(1:nrow(band.year.dt), params6$nr.pix.time.series)
              
              ## Subset dt with band of interest time-series
              band.year.dt.subset <- band.year.dt[pix.to.plot.idx, ]
              band.year.dt.subset <- band.year.dt.subset[order(rank(pixID))]   ## sort pixel to have same order in plots and as a consequence, the same color legend
              dt.band.to.plot <- as.data.table(t(band.year.dt.subset[,-c("pixID", "YrsSince_GrCh", "YrOfChange"), with=FALSE]))   ## remove columns that are not to be melted
              setnames(dt.band.to.plot, band.year.dt.subset[, pixID])  ## final dt has one pixel time-series per column
              dt.band.to.plot[, year:=params6$mapped.years]
              temporal.traj.band.melted <- melt(dt.band.to.plot, id="year")  ## melted dt (repeated entries for each pixel) is the format required by ggplot
              
              ## Time-series of Band 5
              fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("PixTimeSeries_1984_2012_%s_%s_%s_band5.pdf", temporal.data.type, change.type, gsub(" ", "", ecozone.name)), sep='')
              title.str <- sprintf("Band 5 %s time-series:\n%s, %s (%s pix.)", temporal.data.type, change.type, ecozone.name, params6$nr.pix.time.series)
              pdf(fig.name.str)
              theme_set(theme_gray(base_size = 11))
              time.series.band <- ggplot(temporal.traj.band.melted, aes(x=year, y=value, colour=variable, group=variable)) + 
                scale_colour_brewer(palette="Paired") +
                scale_x_continuous(breaks=seq(1984, 2012, by=1), lim=c(1984,2012)) +
                xlab("year") +
                ylab("band 5") +
                coord_cartesian(ylim=params6$ylim.band5) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor.x=element_blank(), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed')) +
                geom_line() +
                ggtitle(title.str)
              print(time.series.band)
              dev.off()
              
            }  ## end if on targ == params6$targ.names[1] & method == params6$methods[1] to plot band 5
            
            
#### MELT FOREST ATTR DATA ---------------------------------------------------
            
            ## Build melted dt to plot time-series of boxplots
            
            ## Get landcover class for every pixel (rows) in each year (columns)
            landcover.dt.subs <- data.table(pixID=yearly.predictions.full[, pixID])  ## initialize dt with pixID
            setkey(landcover.dt.subs, pixID)
            landcover.dt.subs <- merge(landcover.dt.subs, landcover, all=FALSE)  ## merge based on common pixID, as initial "metrics" files have been subset and contain thus less rows
            landcover.matrix <- data.matrix(landcover.dt.subs[, sprintf('%s', params6$mapped.years), with=FALSE])  ## corresponding matrix without pixID
            landcover.prechange.matrix <- landcover.matrix   ## store as a different matrix as landcover.prechange.matrix will be subset when change.type != "always_treed"
            
            if (change.type == "always_treed") {
              
              ## Years since change is a matrix of constant values (50), with some exceptions
              years.since.change.matrix <- matrix(rep(band.year.dt$YrsSince_GrCh, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE) 
              
              xaxis.var <- "sampling_year"  ## dt column name for the variable to plot on the x-axis
              xlabel <- "year"  ## for no change x-axis represents the year 
    
            } else {
              
              ## Retrieve years since change data for every pixel (rows) in each year (columns)
              years.since.change.matrix <- matrix(nrow=nrow(band.year.dt), ncol=length(params6$mapped.years))  ## initialize empty matrix
              for (yr in 1:length(params6$mapped.years)) {
                year <- params6$mapped.years[yr]
                years.since.change.matrix[, yr] <- year - band.year.dt[, YrOfChange]
              }
              
              landcover.prechange.matrix[years.since.change.matrix >= 0] <- 0   ## set to zero landcover values for years after change, including year of change
              
              xaxis.var <- "year_since_ch"
              xlabel <- sprintf("year since %s", change.type)   ## for change x-axis represents the year since greatest change
              
            }
            
            ## Get majority class in pre-change years
            prechange.class <- apply(landcover.prechange.matrix, 1, maj_vote_non_zero)  # apply row-wise the majority vote function that returns most frequent label, excluding zeros 
            prechange.class.matrix <- matrix(rep(prechange.class, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE)  ## replicate vector columnwise to match shape of other matrices
            
            ## Stack together flattened matrices
            dt.yearly.pred.melted <- data.table(sampling_year=c( matrix(rep(params6$mapped.years, each=nrow(yearly.predictions.full)), nrow=nrow(yearly.predictions.full), byrow=FALSE) ),   ## c() reshapes matrix to vector
                                                pix_ID=c( matrix(rep(yearly.predictions.full$pixID, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE) ),
                                                year_since_ch=c( years.since.change.matrix ),    
                                                forest_attr=c( data.matrix(yearly.predictions.full[, sprintf('%s', params6$mapped.years), with=FALSE]) ),
                                                landcover_class=c( landcover.matrix ),  
                                                pre_change_class=c( prechange.class.matrix )
            )
            
            dt.yearly.pred.melted <- dt.yearly.pred.melted[year_since_ch >= params6$start.change.plots]  ## threshold the 2-col dt to keep only attribute values after year 0 (if params6$start.change.plots = 0). Negative values of year_since_ch represent "time before change". In the case of Always_treed this results in keeing all the rows, as time since change values are 50, which should always be > params6$start.change.plots 
            dt.yearly.pred.melted <- dt.yearly.pred.melted[landcover_class %in% params6$forest.classes$AllForest]  ## Subset to plot only forested pixels
            dt.yearly.pred.melted[, year_since_ch:=as.factor(year_since_ch)]  ## set "year_since_ch" as factor to have one boxplot per year
            dt.yearly.pred.melted[, sampling_year:=as.factor(sampling_year)]  ## set "sampling_year" as factor to have one boxplot per year
            
#### BOXPLOTS ---------------------------------------------------
            
            for (i in 1:length(params6$forest.classes)) {
                
              forest.class.name <- names(params6$forest.classes)[i]
              forest.labels <- params6$forest.classes[[i]]
              
              print(sprintf("                     Forest type: %s", forest.class.name))
              
              
              ## Subset to keep only pixels of a certain pre-change class
              dt.yearly.pred.melted.to.plot <- dt.yearly.pred.melted[pre_change_class %in% forest.labels]  
              
              counts.per.year <- data.frame(table(dt.yearly.pred.melted.to.plot[, xaxis.var, with=FALSE]))  ## number of pixels in each yearly bin to plot reference histograms
              nr.plotted.pix <- length(unique(dt.yearly.pred.melted.to.plot$pix_ID))  ## number of pixels that participate to this plot
              
              ## Single code block to produce plots for all change.type options (change and no-change classes)
              fig.name.str <- file.path(FullBoxPlots.subdir, sprintf("Boxplots_1984_2012_%s_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, forest.class.name, gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
              title.str <- sprintf("%s %s boxplots:\n%s, %s, %s, %s (%s pix.)", params6$targ.names.plots[idx.targ], temporal.data.type, change.type, forest.class.name, ecozone.name, method, nr.plotted.pix)
              pdf(fig.name.str)
                theme_set(theme_gray(base_size = 11))
                
                ## Boxplots over time, either absolute time (year) or relative time (since change)
                p1 <- ggplot(dt.yearly.pred.melted.to.plot) +
                      geom_boxplot(aes_string(x=xaxis.var, y="forest_attr"), outlier.shape=NA, notch=F) +
                      ggtitle(title.str) +
                      coord_cartesian(ylim=ylim.targ.abs) +
                      xlab(xlabel) +
                      ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
                      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
              
                ## Histogram with number of pixels in each time since change bin
                p2 <- ggplot(counts.per.year, aes(x=Var1, y=Freq)) +
                      geom_bar(stat="identity") +
                      xlab(xlabel) +
                      ylab("pixel count") +
                      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
                
              if (params6$plot.count.hist) {
                gp1 <- ggplotGrob(p1)
                gp2 <- ggplotGrob(p2)
                grid.draw(rbind(gp1, gp2))  ## to have plots neatly aligned under each other based on y-axis
              } else {
                grid.arrange(p1, nrow=1)
              }
              dev.off()
            
            }
            
            
            #### To check values plotted #### 
            # quantile(yearly.predictions.full[, '2012', with=F][[1]], probs=c(0.1, 0.25, 0.5, 0.75, 0.9))  ## list subsetting of the one-column data.table gives a vector
            # quantile(dt.yearly.pred.melted.to.plot[year_since_ch==-19, forest_attr], probs=c(0.1, 0.25, 0.5, 0.75, 0.9))  ## list subsetting of the one-column data.table gives a vector
            
#### SINGLE-PIXEL TIME-SERIES ---------------------------------------------------
            
            ## Subset dt with predicted resp. var. time-series and transpose to match ggplot time-series format
            dt.final.subset <- yearly.predictions.full[pix.to.plot.idx, ]
            dt.final.subset <- dt.final.subset[order(rank(pixID))]
            setkey(dt.final.subset, pixID)
            dt.final.subset <- merge(dt.final.subset, band.year.dt[, .(pixID, YrOfChange)], all=FALSE)  ## attach YrOfChange...
            dt.final.subset[, pixID:=sprintf("%s_%s", pixID, YrOfChange)]   ## ...and include it in pixID to appear in plot legend
            dt.to.plot <- as.data.table(t(dt.final.subset[,-c("pixID", "YrOfChange"), with=FALSE]))  ## transpose and filter data to meet ggplot format
            setnames(dt.to.plot, dt.final.subset[, pixID])
            dt.to.plot[, year:=params6$mapped.years]
            temporal.traj.melted <- melt(dt.to.plot, id="year")
            
            ## To manually set colors
            # scale.classes <- rep_len(c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"), params6$nr.pix.time.series)
            
            ## Create common ggplot object
            time.series.common <- ggplot(temporal.traj.melted, aes(x=year, y=value, colour=variable, group=variable)) + 
              scale_colour_brewer(palette="Paired") +
              scale_x_continuous(breaks=seq(1984, 2012, by=1), lim=c(1984,2012)) +
              xlab("year") +
              ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
              coord_cartesian(ylim=ylim.targ.abs) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor.x=element_blank(), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
            
            ## Raw pixel time-series
            fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("PixTimeSeries_1984_2012_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
            title.str <- sprintf("%s %s time-series:\n%s, %s, %s (%s pix.)", params6$targ.names.plots[idx.targ], temporal.data.type, change.type, ecozone.name, method, params6$nr.pix.time.series)
            pdf(fig.name.str)
            theme_set(theme_gray(base_size = 11))
              time.series <- time.series.common + 
                              geom_line() + 
                              ggtitle(title.str)
            print(time.series)
            dev.off()
            
            ## Smoothed pixel time-series
            if (params6$plot.smoothed.ts) {
              fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("PixTimeSeries_1984_2012_smoothed_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
              title.str <- sprintf("%s %s smoothed time-series:\n%s, %s, %s (%s pix.)", params6$targ.names.plots[idx.targ], temporal.data.type, change.type, ecozone.name, method, params6$nr.pix.time.series)
              pdf(fig.name.str)
              theme_set(theme_gray(base_size = 11))
                time.series <- time.series.common + 
                                stat_smooth(method = "loess", formula = y ~ x, se = F) + 
                                ggtitle(title.str)
              print(time.series)
              dev.off()
            }
            
            
          } ## end for on params6$methods
        
        } ## end for on params6$targ.names.lg
      
      } ## end for on params6$ecozone.code
      
    }  ## end for on params6a$change.types
    
  }  ## end for on params6$temporal.data.types

}  ## end if on params6$check.temporal.full

## delete temp folder with temp rasters
unlink(raster_tmp_dir, recursive = T, force = T)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog6, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
