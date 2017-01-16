## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog6_temporal_analysis.R                             
## Objective: Plot time-series of mapped attributes with boxplots of an ensemble of pixels sampled over large areas or with lines for individual pixel

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
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

## READS:
# - "<UTMzone>_pt_centerpt_training_validation.shp": (from prog1) point shp with centerpoints of selected polygons
# - "<UTMzone>_poly_training_validation.shp": (from prog1) polygon shp with selected polygons (after recheck for 9-plots/polyg after subsetting wrt availability of both LiDAR and forest attributes)

## WRITES:
# - "poly_training_validation_exvars_extract.csv": (for prog3a and prog3) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

#### INIT --------------------------------------------------------------------

print('Prog5: temporal analysis')

rm(list=ls()) ## clear all variables

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params2.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params3.Rdata"
load(param_file)

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_params4a.Rdata"
load(param_file)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")

## Update results directory to point to the FINAL results one
base_results_dir <- file.path(base_dir, "FINAL_RUN", "results", fsep=.Platform$file.sep)

#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params5 <- list()

params5$subsetting <- F    ## to subset the full map csvs (for debugging in development phase)
params5$nr.pts.subset <- 2000
params5$check.temporal.val <- F  ## wheter to run the analysis over time on val set
params5$temporal.val.data.source <- "CSV"  ## whether to read val set values from CSVs produced by sampling with IDL scripts (option "CSV") or from the raw maps in R using the shapefiles (option "MAP")
params5$temp.val.sampling.shp <- "cpt"    ## sampling unit (to be set as "poly") to read back map values on validation plots, to be set as "cpt" to speed things up in development 
params5$check.temporal.full <- T  ## wheter to run the analysis over time on full maps
params5$temporal.data.types <- c("raw", "fitted")  ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends
# params5$temporal.data.types <- c("fitted")
params5$methods <- c("RF", "YAI")   ## for which methods to run the analysis, can accept both c("RF", "YAI") 
params5$mapped.years <- seq(from=1984, to=2012)   ## only used for val boxplots when reading back the images
# params5$mapped.years <- seq(from=1984, to=1996)   
paramsGL$zones <- "UTM13S"   ## redefine UTM zones to consider for params5$temporal.val.data.source <- "MAP" case
# params5$mapped.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West",
#                                  "Hudson Plains", "Montane Cordillera", "Pacific Maritime", "Taiga Cordillera",
#                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")
# params5$ecozone.codes <- c(2, 3, 4, 5, 6, 7, 9, 11, 15, 16, 17, 18)
params5$mapped.ecozones <- list("Boreal Plains", "Boreal Shield West", "Taiga Shield West")
params5$ecozone.codes <- c(4, 6, 18)
params5$targ.names <- c("elev_p95", "percentage_first_returns_above_2m")
params5$targ.units <- c("m", "%")
params5$targ.names.plots <- c("elev_p95", "cover_2m")
params5$ylim.targ$absolute$elev_p95 <- c(0, 28)
params5$ylim.targ$absolute$cover_2m <- c(0, 100)
params5$ylim.targ$difference$elev_p95 <- c(-11, 11)
params5$ylim.targ$difference$cover_2m <- c(-40, 40)
params5$ylim.band5 <- c(400, 1800)

params5$nr.pix.time.series <- 12

param_file_prog5 = file.path(base_wkg_dir, 'AllUTMzones_params5.Rdata', fsep = .Platform$file.sep)
save(params5, file = param_file_prog5)

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

## Subdirectory to save model assessment plots
Temporal.subdir <- file.path(base_figures_dir, "TemporalAnalyses", fsep = .Platform$file.sep)
if (! file.exists(Temporal.subdir)){dir.create(Temporal.subdir, showWarnings = F, recursive = T)}

#### 1) BOXPLOTS ON VAL SET ---------------------------------------------------

## Run only if we check for temporal boxplots VS 2010 observed values
if (params5$check.temporal.val) {  
  
  print('Validation set boxplots:')
  
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
  
  if (params5$temporal.val.data.source == "MAP") {
  
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
      for (yr in params5$mapped.years) {
        cmd <- sprintf('map%s <- raster(file.path("%s", "UTM_results", "RF_corrected", "UTM_%s_elev_p95_RF_%s.dat", fsep = .Platform$file.sep))', yr, Landsat_dir, zone.nr, yr)  
        eval(parse(text=cmd))
      }
      
      vars.years <- c(params5$mapped.years)  
      str.years.vect <- sprintf("map%s", params5$mapped.years)
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
        
        ## Selection of the sampling unit
        if (params5$temp.val.sampling.shp == 'poly') {   
          sampling.shp <- poly.validation[indPart, ]   ## for actual study we use polygons
        } else if (params5$temp.val.sampling.shp == 'cpt') { 
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
    
    setnames(yearly.predictions.R, old=sprintf("X%s", params5$mapped.years), new=as.character(params5$mapped.years))
    yearly.predictions.R <- data.table(FCID=yearly.predictions.R[, FCID], yearly.predictions.R[, as.character(params5$mapped.years), with=FALSE]/1000)  ## convert back to proper values
    
  }  ## end if params5$temporal.val.data.source == "MAP"
  
  
  
#### PLOT BOXPLOTS ---------------------------------------------------------------------
  
  for (targ in params5$targ.names) {
    
    print(sprintf("  Targ: %s", targ))
    
    idx.targ <- which(targ == params5$targ.names)   ## numeric index useful to select elements of string vectors
    
    ## set of target-specific parameters
    if (targ == "elev_p95") {
      observed.val.col <- "O_el_p95"  ## name of the column with the observed values in plot.based.predictions@data
      ylim.targ.abs <- params5$ylim.targ$absolute$elev_p95
      ylim.targ.diff <- params5$ylim.targ$difference$elev_p95
    } else if (targ == "percentage_first_returns_above_2m") {
      observed.val.col <- "O_p_1r_2m"
      ylim.targ.abs <- params5$ylim.targ$absolute$cover_2m
      ylim.targ.diff <- params5$ylim.targ$difference$cover_2m
    }
    
    for (method in params5$methods) { 
      
      print(sprintf("    Method: %s", method))
      
      print("      Ch_attr: No change")
      
      if (params5$temporal.val.data.source == "CSV") {
        file.name <- file.path(val.time.series.folder.name, sprintf("time_series_val_set_%s_%s.csv", method, targ), fsep = .Platform$file.sep)
        yearly.predictions.CSV <- fread(file.name, header = T)
        yearly.predictions.CSV[, FCID:=plots.info$FCID[plots.info$TV=="VALIDATION"]]
        yearly.predictions <- yearly.predictions.CSV
      } else if (params5$temporal.val.data.source == "MAP") {
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
        fig.name.str <- file.path(base_figures_dir, "MapCheck_CAN_level", "Plot_VS_Pixel_CAN", sprintf("Plot_VS_Pixel_CAN_no_weigh_avg_%s_%s.pdf", params5$targ.names.plots[idx.targ], method), sep='')
        pdf(fig.name.str)
           p <- plot_colorByDensity(merged.pred$PlotPred, merged.pred$PixelPred, xlab="plot prediction", ylab="pixel prediction", main=sprintf("%s [%s]: R^2=%.3f on %s validation points", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ], R2coeff, nrow(merged.pred))) +
                geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line                        
        print(p)
        dev.off()
        
        ## On UTM zones 12S and 13S only (to compare with reading back preliminary maps)
        subset.UTM12S13S <- X.trn.val.raw.w.ECO %>%
          filter(TV=="VALIDATION" & (UTMzone=="13S" | UTMzone=="12S") )
        subset.UTM12S13S <- as.data.table(subset.UTM12S13S)
        merged.pred.subs <- merge(subset.UTM12S13S, merged.pred, by="FCID")
        R2coeff <- cor(merged.pred.subs$PlotPred, merged.pred.subs$PixelPred)**2
        fig.name.str <- file.path(base_figures_dir, "MapCheck_CAN_level", "Plot_VS_Pixel_UTM12S13S", sprintf("Plot_VS_Pixel_UTM12S13S_no_weigh_avg_%s_%s.pdf", params5$targ.names.plots[idx.targ], method), sep='')
        pdf(fig.name.str)
           p <- plot_colorByDensity(merged.pred.subs$PlotPred, merged.pred.subs$PixelPred, xlab="plot prediction", ylab="pixel prediction", main=sprintf("%s [%s]: R^2=%.3f on %s validation points", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ], R2coeff, nrow(merged.pred.subs))) +
                geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line
        print(p)
        dev.off()
        
      }
      
      ## merge the three dataframes by FCID to have corresponding values of pixel and plot based predictions
      X.yearly.predictions <- merge(X.trn.val.raw.w.ECO, yearly.predictions, by="FCID")
      dt.with.year <- merge(as.data.table(plot.based.predictions@data[, c("FCID", observed.val.col)]), X.yearly.predictions, by="FCID")
      
      ## Predicted boxplots over time on validation set + observed boxplot in 2010
      for (ez in 1:length(params5$mapped.ecozones)) {
        
        ecozone.name <- params5$mapped.ecozones[ez]
        
        ## Skip if the ecozone in not among those found in the plots 
        if ( !ecozone.name %in% unique(dt.with.year$ECOZONE) ) { next } 
        
        ## Select only no-change plots of the ecozone of interest
        dt.final <- dt.with.year %>%
                    filter(ECOZONE==ecozone.name, YrsSince_GrCh==50, Ch_attr==0)
        dt.final <- as.data.table(dt.final)
        dt.final <- dt.final[, c(params5$mapped.years, observed.val.col), with=FALSE]
        
        ## Rename columns to reflect predicted (P_) and observed (O_) values
        colnames(dt.final) <- c( sprintf("P_%s", colnames(dt.final)[1:length(colnames(dt.final))-1]), "O_2010")
        
        ## Computing residuals wrt to 2010 observed value
        dt.final.point.diff <- as.data.table(apply(dt.final, 2, FUN=function (x) x-dt.final$O_2010 ))
        
        yearly.predictions <- melt(dt.final)
        yearly.predictions.point.diff <- melt(dt.final.point.diff)
        
        ## Boxplot of predicted values compared to observed values in 2010 
        fig.name.str <- file.path(ValBoxPlots.subdir, sprintf("NoChange_ValSet_%s_%s_%s_1984_2012.pdf", gsub(" ", "", ecozone.name), params5$targ.names.plots[idx.targ], method), sep='')
        title.str <- sprintf("Val. set:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, nrow(dt.final))
        pdf(fig.name.str)
        theme_set(theme_gray(base_size = 14))
        p <- ggplot(yearly.predictions) +
          geom_boxplot(aes(x=variable, y=value), outlier.shape=NA, notch=F) +
          ggtitle(title.str) +
          coord_cartesian(ylim=ylim.targ.abs) +   ## don't use ylim() as it will change shape bc it remove points beyond limit !!!
          xlab("year") +
          ylab(sprintf("%s [%s]", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ])) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
        print(p)
        dev.off()
        
        # Difference with observed 2010 value (residuals wrt 2010)
        fig.name.str <- file.path(ValBoxPlots.subdir, sprintf("NoChange_ValSet_Residuals_%s_%s_%s_1984_2012.pdf", gsub(" ", "", ecozone.name), params5$targ.names.plots[idx.targ], method), sep='')
        title.str <- sprintf("Per-sample difference (residual) wrt 2010 \n Val. set:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, nrow(dt.final))
        pdf(fig.name.str)
        theme_set(theme_gray(base_size = 14))
        p <- ggplot(yearly.predictions.point.diff) +
          geom_boxplot(aes(x=variable, y=value), outlier.shape=NA, notch=F) +
          ggtitle(title.str) +
          coord_cartesian(ylim=ylim.targ.diff) +   ## don't use ylim() as it will change shape bc it remove points beyond limit !!!
          xlab("year") +
          ylab(sprintf("residual %s [%s]", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ])) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
        print(p)
        dev.off()
          
      }   ## end for on params5$ecozone.code
      
    }  ## end for on params5$methods
     
  }   ## end for on params5$targ.names
  
}  ## end if params5$check.temporal.val

#### 2) PLOTS ON FULL MAP ---------------------------------------------------
  
if (params5$check.temporal.full) {
  
  for (temporal.data.type in params5$temporal.data.types) {
  
    print(sprintf("Temporal data: %s", temporal.data.type))
    
    print('Full map plots:')
    
    if (temporal.data.type == "fitted") {
      samples.dir <- file.path(params4a$temporal.base.dir, 'samples_fitted', fsep = .Platform$file.sep) 
      predicted.values.dir <- 'predicted_values_fitted'
      predictors.dir <- file.path(params4a$temporal.base.dir, 'metrics_fitted', fsep = .Platform$file.sep) 
      FullBoxPlots.subdir <- file.path(Temporal.subdir, "FullBoxPlots_fitted", fsep = .Platform$file.sep)
      FullPixelTimeSeries.subdir <- file.path(Temporal.subdir, "FullPixelTimeSeries_fitted", fsep = .Platform$file.sep)
    } else if (temporal.data.type == "raw") {
      samples.dir <- file.path(params4a$temporal.base.dir, 'samples', fsep = .Platform$file.sep) 
      predicted.values.dir <- 'predicted_values'
      predictors.dir <- file.path(params4a$temporal.base.dir, 'metrics', fsep = .Platform$file.sep) 
      FullBoxPlots.subdir <- file.path(Temporal.subdir, "FullBoxPlots", fsep = .Platform$file.sep)
      FullPixelTimeSeries.subdir <- file.path(Temporal.subdir, "FullPixelTimeSeries", fsep = .Platform$file.sep)
    }
    if (! file.exists(FullBoxPlots.subdir)){dir.create(FullBoxPlots.subdir, showWarnings = F, recursive = T)}
    if (! file.exists(FullPixelTimeSeries.subdir)){dir.create(FullPixelTimeSeries.subdir, showWarnings = F, recursive = T)}
    folder.name <- file.path(Landsat_dir, "UTM_results", "random_samples_always_forested", predicted.values.dir, fsep = .Platform$file.sep)  
  
    for (ez in 1:length(params5$mapped.ecozones)) {
      
      ecozone.name <- params5$mapped.ecozones[ez]
      ecozone.code <- params5$ecozone.codes[ez]
      
      print(ecozone.name)
      
      for (targ in params5$targ.names) {
          
        print(sprintf("  Targ: %s", targ))
        
        idx.targ <- which(targ == params5$targ.names)   ## numeric index useful to select elements of string vectors
  
        for (method in params5$methods) { 
          
          print(sprintf("    Method: %s", method))
  
  
#### NO CHANGE PIXELS -----------------------------------------------------------
          
          print("      Ch_attr: No change")
          
          ## Yearly boxplots in persistent no change areas, full map
          cmd <- sprintf('file.name <- file.path(folder.name, "persistant_treed_%s_%s_%s.csv", fsep = .Platform$file.sep)', method, targ, gsub(' ', '', ecozone.name))  
          eval(parse(text=cmd))
          if (!file.exists(file.name)){ next } ## if file does not exist skip to next round of loop on targ
   
          ## set of target-specific parameters
          if (targ == "elev_p95") {
            dt.targ.col.name <- "p95"  ## only used for change csvs
            ylim.targ.abs <- params5$ylim.targ$absolute$elev_p95
            alpha.targ <- 0.2
          } else if (targ == "percentage_first_returns_above_2m") {
            dt.targ.col.name <- "cc"
            ylim.targ.abs <- params5$ylim.targ$absolute$cover_2m
            alpha.targ <- 0
          }
          
          ## read csvs (and subset them if specified)
          if (!params5$subsetting) {  ## if not subsetting read full file...
            yearly.predictions.full <- fread(file.name, header = T)
          } else {  ## ...else specify a subset directory (_subset)
            folder.name.subsets <- paste0(folder.name, '_subset')
            if (! file.exists(folder.name.subsets)){dir.create(folder.name.subsets, showWarnings = F, recursive = T)}
            cmd <- sprintf('file.name.subs <- file.path(folder.name.subsets, "persistant_treed_%s_%s_%s.csv", fsep = .Platform$file.sep)', method, targ, gsub(' ', '', ecozone.name))  
            eval(parse(text=cmd))
            if (file.exists(file.name.subs)) {  ## that if it already exists...
              yearly.predictions.full <- fread(file.name.subs, header = T)   ## ...we read
            } else {  ## ...otherwise we read the full one, subset it and write it 
              ## read file into a data.table (faster than readCSV)
              yearly.predictions.full <- fread(file.name, header = T)
              colnames.order <- colnames(yearly.predictions.full)   ## save initial column names as the order is changed by sample()
              set.seed(paramsGL$global.seed)
              yearly.predictions.full <- yearly.predictions.full[,.SD[sample(.N, params5$nr.pts.subset)], by=ecozone]
              setcolorder(yearly.predictions.full, colnames.order)  ## to reput ecozone at the end
              fwrite(yearly.predictions.full, file.name.subs)
            }
          }
          
          ## Band 5 data time-series for comparison with time-series of predicted data
          if (targ == params5$targ.names[1] & method == params5$methods[1]) {  ## to do this plot only once per ecozone
            
            file.name <- file.path(samples.dir, sprintf("%s_%s.csv", ecozone.code, ecozone.name), fsep = .Platform$file.sep)
            pix.IDs <- fread(file.name, header = F)
            
            band.year.dt <- data.table(pixID=yearly.predictions.full$pixID, idx=seq(from=1, to=nrow(yearly.predictions.full)))
            setkey(band.year.dt, pixID)
            for (yr in 1:length(params5$mapped.years)) {
              
              ## Load dataset from CSV
              file.name <- file.path(predictors.dir, sprintf("%s_%s_%s.csv", ecozone.code, ecozone.name, params4a$mapped.years[yr]), fsep = .Platform$file.sep)
              yearly.predictors <- fread(file.name, header = T)
              yearly.predictors[, pixID:=sprintf('%s_%d', pix.IDs[[1]], pix.IDs[[2]])]
              yearly.predictors[, index:=NULL]
              setkey(yearly.predictors, pixID) 
              
              band.year.dt <- merge(band.year.dt, yearly.predictors[, .(pixID, b5)], all=FALSE)
              colnames(band.year.dt)[ncol(band.year.dt)] <- as.character(params5$mapped.years[yr])
              
            }
            
            setkey(band.year.dt, idx)  ## reset order
            band.year.dt[, "idx":=NULL]   ## ...and then remove it

            ## Get pixel indices to plot in time-series
            set.seed(paramsGL$global.seed)
            pix.to.plot.idx <- sample(1:nrow(band.year.dt), params5$nr.pix.time.series)
            
            ## Subset dt with band of interest time-series
            band.year.dt.subset <- band.year.dt[pix.to.plot.idx, ]
            dt.band.to.plot <- as.data.table(t(band.year.dt.subset[,-"pixID", with=FALSE]))
            setnames(dt.band.to.plot, band.year.dt.subset[, pixID])
            dt.band.to.plot[, year:=params5$mapped.years]
            temporal.traj.band.melted <- melt(dt.band.to.plot, id="year")
            
            ## Band of interest time-series
            fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("NoChange_PixTimeSeries_Full_%s_band5_1984_2012.pdf", gsub(" ", "", ecozone.name)), sep='')
            title.str <- sprintf("Band 5 time-series:\n%s, perm. no change (%s pix.)", ecozone.name, params5$nr.pix.time.series)
            pdf(fig.name.str)
            theme_set(theme_gray(base_size = 14))
            time.series.band <- ggplot(temporal.traj.band.melted, aes(x=year, y=value, colour=variable, group=variable)) + 
              scale_colour_brewer(palette="Paired") +
              scale_x_continuous(breaks=seq(1984, 2012, by=1), lim=c(1984,2012)) +
              xlab("year") +
              ylab("band 5") +
              coord_cartesian(ylim=params5$ylim.band5) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor.x=element_blank(), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed')) +
              geom_line() + 
              ggtitle(title.str)
            print(time.series.band)
            dev.off()
            
          }  ## end if on targ == params5$targ.names[1] & method == params5$methods[1] to plot band 5
          
          
          dt.final.full <- yearly.predictions.full[, -"ecozone", with=FALSE]
          dt.yearly.pred.by.eco.full <- melt(dt.final.full)
          
          fig.name.str <- file.path(FullBoxPlots.subdir, sprintf("NoChange_Full_%s_%s_%s_1984_2012.pdf", 
                                                                       gsub(" ", "", ecozone.name), 
                                                                       params5$targ.names.plots[idx.targ], 
                                                                       method), sep='')
          title.str <- sprintf("Maps:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, nrow(dt.final.full))
          pdf(fig.name.str)
            theme_set(theme_gray(base_size = 14))
            p <- ggplot(dt.yearly.pred.by.eco.full) +
              geom_boxplot(aes(x=variable, y=value), outlier.shape=NA, notch=F) +
              # scale_x_discrete(breaks=seq(1984, 2012, by=2), lim=c(1984,2012)) +
              ggtitle(title.str) +
              coord_cartesian(ylim=ylim.targ.abs) +
              xlab("year") +
              ylab(sprintf("%s [%s]", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ])) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
          print(p)
          dev.off()
          
          ## to check values plotted
          # quantile(yearly.predictions.full[, '1984', with=F][[1]], probs=c(0.1, 0.25, 0.5, 0.75, 0.9))  ## list subsetting of the one-column data.table gives a vector
  
          ## Subset dt with predicted resp. var. time-series and transpose to match ggplot time-series format
          dt.final.subset <- dt.final.full[pix.to.plot.idx, ]
          dt.to.plot <- as.data.table(t(dt.final.subset[,-"pixID", with=FALSE]))
          setnames(dt.to.plot, dt.final.subset[, pixID])
          dt.to.plot[, year:=params5$mapped.years]
          temporal.traj.melted <- melt(dt.to.plot, id="year")
          
          ## To manually set colors
          # scale.classes <- rep_len(c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"), params5$nr.pix.time.series)
          
          ## Create common ggplot object
          time.series.common <- ggplot(temporal.traj.melted, aes(x=year, y=value, colour=variable, group=variable)) + 
            scale_colour_brewer(palette="Paired") +
            scale_x_continuous(breaks=seq(1984, 2012, by=1), lim=c(1984,2012)) +
            xlab("year") +
            ylab(sprintf("%s [%s]", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ])) +
            coord_cartesian(ylim=ylim.targ.abs) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor.x=element_blank(), panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed'))
          
          ## Raw pixel time-series
          fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("NoChange_PixTimeSeries_Full_%s_%s_%s_1984_2012.pdf", 
                                                                 gsub(" ", "", ecozone.name), 
                                                                 params5$targ.names.plots[idx.targ], 
                                                                 method), sep='')
          title.str <- sprintf("Pixel time-series:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, params5$nr.pix.time.series )
          pdf(fig.name.str)
          theme_set(theme_gray(base_size = 14))
            time.series <- time.series.common + 
                            geom_line() + 
                            ggtitle(title.str)
          print(time.series)
          dev.off()
          
          ## Smoothed pixel time-series
          fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("NoChange_PixTimeSeries_Smoothed_Full_%s_%s_%s_1984_2012.pdf", 
                                                                        gsub(" ", "", ecozone.name), 
                                                                        params5$targ.names.plots[idx.targ], 
                                                                        method), sep='')
          title.str <- sprintf("Smoothed pixel time-series:\n%s, %s, perm. no change (%s pix.)", ecozone.name, method, params5$nr.pix.time.series )
          pdf(fig.name.str)
          theme_set(theme_gray(base_size = 14))
            time.series <- time.series.common + 
                            stat_smooth(method = "loess", formula = y ~ x, se = F) + 
                            ggtitle(title.str)
          print(time.series)
          dev.off()
          
          
  #### CHANGE PIXELS -----------------------------------------------------------
          
          ## Yearly boxplots in changed areas, full map
            
          ## TO UNCOMMENT AND ADJUST WHEN CHANGE DATA WILL BE AVAILABLE
          # for (chattr in params3$Ch_attr.classes ) {
          #   
          #   cmd <- sprintf('file.name <- file.path("%s", "UTM_results", "time_results", "%s", "%s_postchange_%s_%sk_%s.csv", fsep = .Platform$file.sep)', Landsat_dir, params5$csv.version, gsub(" ", "", chattr), targ, params5$nr.sampled.pix, params5$csv.version)  
          #   eval(parse(text=cmd))
          #   if (!file.exists(file.name)) { next }   ## if file does not exist skip this round of the loop
          #   
          #   print(sprintf("    Ch_attr: %s", chattr))
          #   
          #   yearly.predictions.full <- fread(file.name, header = T)
          #   
          #   dt.yearly.pred.by.eco.full <- yearly.predictions.full[ecozone==params5$ecozone.code[ez], c(dt.targ.col.name, "time_since"), with=F]  ## this format is already "melted"
          #   dt.yearly.pred.by.eco.full <- dt.yearly.pred.by.eco.full[, time_since:=as.factor(time_since)]
          #   df.counts <- data.frame(table(dt.yearly.pred.by.eco.full$time_since))
          #   
          #   ## detect biggest drop in median value and set this as the boundary between mortality and regrowth
          #   cmd <- sprintf("median.over.time <- arrange(dt.yearly.pred.by.eco.full[, median(%s), by=time_since], time_since)", dt.targ.col.name)
          #   eval(parse(text=cmd))
          #   
          #   biggest.drop <- 0
          #   idx.drop <- NA
          #   for (i in seq(from=1, to=nrow(median.over.time)-1, by=1) ) {
          #     drop <- median.over.time[i, 2] - median.over.time[i+1, 2]  ## 2nd column is always the median value, so it is hard-coded
          #     if (drop > biggest.drop) {
          #       biggest.drop <- drop    
          #       idx.drop <- i
          #     }
          #   }
          #   fig.name.str <- file.path(FullBoxPlots.subdir, sprintf("%s_Full_%s_%s_1984_2012_%sKsamples_%s.pdf", gsub(" ", "", chattr), gsub(" ", "", params5$ecozone.names[ez]), params5$targ.names.plots[idx.targ], params5$nr.sampled.pix), params5$csv.version, sep='')
          #   title.str <- sprintf("%s, %s", params5$ecozone.names[ez], chattr)
          #   pdf(fig.name.str)
          #     theme_set(theme_gray(base_size = 14))
          #     
          #     ## boxplots over time since change
          #     p1 <- ggplot(dt.yearly.pred.by.eco.full) +
          #       geom_boxplot(aes_string(x="time_since", y=dt.targ.col.name), outlier.shape=NA, notch=F) +
          #       # scale_x_discrete(breaks=seq(1984, 2012, by=2), lim=c(1984,2012)) +
          #       ggtitle(title.str) +
          #       coord_cartesian(ylim=ylim.targ) +
          #       xlab(sprintf("years since %s", chattr)) +
          #       ylab(sprintf("%s [%s]", params5$targ.names.plots[idx.targ], params5$targ.units[idx.targ])) +
          #       theme(panel.grid.major.y=element_line(colour='gray', linetype = 'dashed'), panel.grid.minor.y=element_line(colour='gray', linetype = 'dashed')) +
          #       geom_rect(data=data.frame(xstart=0, xend=idx.drop+2.5), aes(xmin=xstart, xmax=xend, ymin=0, ymax=ylim.targ[2]), alpha=alpha.targ, fill="red") +  ## box showing mortality
          #       geom_rect(data=data.frame(xstart=idx.drop-0.5, xend=28), aes(xmin=xstart, xmax=xend, ymin=0, ymax=ylim.targ[2]), alpha=alpha.targ, fill="green")   ## box showing growth
          #     
          #     ## number of pixels in each time since change bin
          #     p2 <- ggplot(df.counts, aes(x=Var1, y=Freq)) + 
          #       geom_bar(stat="identity") + 
          #       xlab(sprintf("years since %s", chattr)) +
          #       ylab("pixel count")
          #     
          #   grid.arrange(p1, p2, nrow=2)
          #   dev.off()
          #     
          # } ## end for on params5$Ch_attr.classes
          
        } ## end for on params5$methods
      
      } ## end for on params5$targ.names.lg
    
    } ## end for on params5$ecozone.code
    
  }  ## end for on params5$temporal.data.types

}  ## end if on params5$check.temporal.full

## delete temp folder with temp rasters
unlink(raster_tmp_dir, recursive = T, force = T)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog5, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
