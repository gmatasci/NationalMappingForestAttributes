## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca)       
## File Name: prog6_temporal_analysis.R                             
## Objective: Plot time-series of mapped attributes with boxplots of an ensemble of pixels sampled over large areas or with lines for individual pixel

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# -

## SOLVED:
# -V change back and remove _v2 -- now we have separate folders for each version of the csvs
# -V save separate correlation trends plots -- plots are now produced
# -V Check for 0,0,0,0 in bands values: remove them! -- removed in prog4c with line  yearly.predictors <- yearly.predictors[b1!=0 & b2!=0 & b3!=0 & b4!=0 & b5!=0 & b6!=0]
# -V fix pixID issue band5 wrt rest of plots -- order of sampled pixels is preserved by using merge to keep matching pixel IDs instead of cbinding data.tables which can lead to mixing different orders
# -V implement and check script for new boxplots v3 -- v3 boxplots are more stable bc of the further filtering of the sampled pixels


#### READS/WRITES ------------------------------------------------------------


#### INIT --------------------------------------------------------------------

print('Prog6: temporal analysis')

rm(list=ls()) ## clear all variables

param_file = "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

base_dir <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR"
base_wkg_dir <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg"
base_figures_dir <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/figures"
base_results_dir <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/results"


param_file = "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg/AllUTMzones_params2.Rdata"
load(param_file)

param_file = "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg/AllUTMzones_params3.Rdata"
load(param_file)

param_file = "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg/AllUTMzones_params6a.Rdata"
load(param_file)


source("D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/code/Functions_NatMapping_R.R")


#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params6 <- list()

params6$experiment.name <- "fitted_NO_CHATTR"  ## matching to folder names in E:\NTEMS\UTM_results\2017_02_27_SAMPLES_treed_changes_plots\<ch_attr> and in D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR

# params6$plot.description <- "NTL10_RL3_CANplotStartAt10_1984_2016"
# params6$plot.description <- "LorHt_AgBiom"
params6$plot.description <- "Elevp95_cover2m"

params6$plot.boxplot.counts <- F  ## whether or not to plot boxplot counts

params6$temporal.data.types <- c("fitted")  ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends

params6a$samples.folder.name <- 'samples_u'   ## folder with sample index and UTMzone
params6a$landcover.folder.name <- 'landcover_s'  ## VLCE predicted class to mask and filter time-series

params6$change.types <- c("always_treed", "fire", "harvesting")   ## include all of them or code will raise errors when mapping the 3 trends on same plot
# params6$change.types <- c("fire", "harvesting")   ## include all of them or code will raise errors when mapping the 3 trends on same plot

params6$methods <- c("YAI")   ## for which methods to run the analysis, can accept both c("RF", "YAI") 

params6$mapped.years <- seq(from=1984, to=2016)     ## years to map

params6$forest.type.period <- "post-change"  ## define forest type based on either pre- ("pre-change") or post-change ("post-change") majority class

params6$no.tree.lag <- 10    ## number of years after disturbance during which no trees are expected: in this period we fit 2nd degree polynimial
params6$residual.lag <- 3   ## number of years after disturbance during which there should be no residual trees for the pixel to be classified as regrowing (otherwise, if there are trees in this period we classify as residual)

# params6$forest.classes <- list("AllForest"=c(81, 210, 220, 230), "WetlandTreed"=81, "ConiferousForest"=210, "DeciduousForest"=220, "MixedForest"=230)
params6$forest.classes <- list("AllForest"=c(81, 210, 220, 230), "ConiferousForest"=210, "DeciduousForest"=220)

params6$trend.types <- list("AlwTreeInChange"="alwTree", "Mixed"=c("regrow", "residual"), "Regrow"="regrow", "Residual"="residual")

params6$boxplot.title <- 'short'  ## set to 'long' for testing phase to have a telling text in the figure, set to 'short' when producing figures for paper
# params6$boxplot.title <- 'long'

# params6$start.change.plots <- -50  ## year since gr change at which to start plotting the boxplot time-series (x-axis lower limit)
params6$start.change.plots <- -3
params6$start.regrowth.ecoz.plots <- 0   ## year since gr change at which to start plotting the regrowing phase in the ecozone boxplots time-series
params6$start.regrowth.CAN.plots <- 10  ## year since gr change at which to start plotting the regrowing phase in the Canada-wide (combined ecozones) boxplots time-series

params6$start.val$elev_p95 <- 2   ## start value for elev_p95 modeling data when interpolating first 10 years
params6$start.val$cover_2m <- 5   ## start value for cover_2m modeling data when interpolating first 10 years
params6$start.val$stem_volume <- 2     
params6$start.val$ag_biomass <- 2  
params6$start.val$loreys_height <- 3.5
params6$start.val$basal_area <- 0.8   

params6$boxplots$color$always.treed <- "black"   ## colors for boxplots and lines in Canada-wide plot
params6$boxplots$color$harvesting <- "darkgreen"
params6$boxplots$color$fire <- "red4"
params6$boxplots$alpha <- 0.2    ## shaded area transparency in plots

params6$plot.smoothed.ts <- F    ## whether to produce Loess smoothed plots of single pixel time-series

## List of ecozones to plot: params6$fire.ecozones and params6$harvest.ecozones have to be contained in params6$mapped.ecozones, otherwise errors are raised
params6$mapped.ecozones <- list("Atlantic Maritime"=2, "Boreal Cordillera"=3, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Hudson Plains"=7, "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Cordillera"=15,
                                 "Taiga Plains"=16, "Taiga Shield East"=17, "Taiga Shield West"=18)
params6$harvest.ecozones <- list("Atlantic Maritime"=2, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Plains"=16)

params6$fire.ecozones <- list("Boreal Cordillera"=3, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Hudson Plains"=7, "Montane Cordillera"=9, "Taiga Cordillera"=15,
                                 "Taiga Plains"=16, "Taiga Shield East"=17, "Taiga Shield West"=18)
 
# params6$mapped.ecozones <- list("Atlantic Maritime"=2, "Montane Cordillera"=9, "Pacific Maritime"=11)
# 
# params6$harvest.ecozones <- list("Atlantic Maritime"=2, "Pacific Maritime"=11)
# 
# params6$fire.ecozones <- list("Montane Cordillera"=9)

## Parameters of plot axis
params6$targ.names <- c("elev_p95", "percentage_first_returns_above_2m")
params6$targ.units <- c("m", "%")
params6$targ.names.plots <- c("elev_p95", "cover_2m")
params6$ylim.targ$absolute$elev_p95$NONBOREAL <- c(0, 40)
params6$ylim.targ$absolute$elev_p95$BOREAL <- c(0, 28)
params6$ylim.targ$absolute$percentage_first_returns_above_2m$NONBOREAL <- c(0, 100)
params6$ylim.targ$absolute$percentage_first_returns_above_2m$BOREAL <- c(0, 100)

# params6$targ.names <- c("loreys_height", "total_biomass")
# params6$targ.units <- c("m", "t/ha")
# params6$targ.names.plots <- c("loreys_height", "ag_biomass")
# params6$ylim.targ$absolute$loreys_height$NONBOREAL <- c(0, 35)
# params6$ylim.targ$absolute$loreys_height$BOREAL <- c(0, 22)
# params6$ylim.targ$absolute$total_biomass$NONBOREAL <- c(0, 500)
# params6$ylim.targ$absolute$total_biomass$BOREAL <- c(0, 200)

# params6$targ.names <- c("basal_area", "gross_stem_volume")
# params6$targ.units <- c("m^2/ha", "m^3/ha")
# params6$targ.names.plots <- c("basal_area", "stem_volume")
# params6$ylim.targ$absolute$basal_area$NONBOREAL <- c(0, 90)
# params6$ylim.targ$absolute$basal_area$BOREAL <- c(0, 60)
# params6$ylim.targ$absolute$gross_stem_volume$NONBOREAL <- c(0, 1800)
# params6$ylim.targ$absolute$gross_stem_volume$BOREAL <- c(0, 800)

params6$ylim.band5 <- c(200, 3000)

params6$nr.pix.time.series <- 12    ## number of pixels to consider in the single-pixel time-series

param_file_prog6 = file.path(base_wkg_dir, 'AllUTMzones_params6.Rdata', fsep = .Platform$file.sep)
save(params6, file = param_file_prog6)

#### LOAD PACKAGES ----------------------------------------------------------

## install data.table version with fwrite
# install.packages("data.table", repos = "https://Rdatatable.github.io/data.table", type = "source")

list.of.packages <- c("gridExtra",
                      # "dplyr",    ## to be loaded before foreach to avoid "assertion failed" errors"ggplot2",
                      "grid",
                      "reshape2",
                      "rgdal",
                      "raster",
                      "GSIF",
                      "sp",
                      "plyr", 
                      "landsat",
                      "lattice",
                      "lubridate", 
                      "doParallel",
                      "foreach",
                      "data.table",
                      "methods",
                      "egg",
                      "trend"
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


#### PLOTS ON FULL MAP ---------------------------------------------------


print(sprintf('Full map plots: %s', params6$experiment.name))

## Loop on type of temporal data: raw or fitted
for (temporal.data.type in params6$temporal.data.types) {
  
  print(sprintf("Temporal data: %s", temporal.data.type))
  
  list.dt.CAN.plot <- list()  ## initialize here list containing dt for CANADA plot fire vs harvest
  list.trend.type.counts <- list()  ## initialize here list containing stats on how many pixels per trend type do we have
  
  ## Initialize dynamic list of empty tables to save number of plotted pixels
  list.nr.plotted.pix <- list()
  for (i in 1:length(params6$forest.classes)) { 
    forest.class.name <- names(params6$forest.classes)[i]
    for (t in 1:length(params6$trend.types)) {
      trend.type.name <- names(params6$trend.types)[t]
      empty.dt <- as.data.table( data.frame(matrix(nrow=max(unlist(params6$mapped.ecozones)), ncol=length(params6$change.types))) )
      colnames(empty.dt) <- params6$change.types
      empty.dt <- empty.dt[, lapply(.SD, as.numeric)]
      list.nr.plotted.pix[[sprintf('%s_%s', forest.class.name, trend.type.name)]] <- empty.dt
    }
  }
  
  ## Loop on types of change 
  for (change.type in params6$change.types) {
    
    print(sprintf("    Ch_attr: %s", change.type))
    
    samples.dir <- file.path(params6a$temporal.base.dir, change.type, params6a$samples.folder.name, fsep = .Platform$file.sep)
    landcover.dir <- file.path(params6a$temporal.base.dir, change.type, params6a$landcover.folder.name, fsep = .Platform$file.sep)
    
    ## Choose type of temporal data
    if (temporal.data.type == "fitted") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics_fitted', fsep = .Platform$file.sep)
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, sprintf('predicted_values_%s', params6$experiment.name), fsep = .Platform$file.sep)
      Temporal.subdir <- file.path(sprintf('%s_%s', base_figures_dir, params6$experiment.name), sprintf('TmpAn_%s', params6$plot.description), fsep = .Platform$file.sep)
    } else if (temporal.data.type == "raw") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics', fsep = .Platform$file.sep) 
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, 'predicted_values', fsep = .Platform$file.sep) 
      Temporal.subdir <- file.path(base_figures_dir, sprintf('TmpAn_%s', params6$plot.description), fsep = .Platform$file.sep)
    }
    FullBoxPlots.subdir <- file.path(Temporal.subdir, "FullBoxPlots", change.type, fsep = .Platform$file.sep)
    FullPixelTimeSeries.subdir <- file.path(Temporal.subdir, "FullPixelTimeSeries", change.type, fsep = .Platform$file.sep)
    if (! file.exists(FullBoxPlots.subdir)){dir.create(FullBoxPlots.subdir, showWarnings = F, recursive = T)}
    if (! file.exists(FullPixelTimeSeries.subdir)){dir.create(FullPixelTimeSeries.subdir, showWarnings = F, recursive = T)}
    
    ## If source folder does not exist, skip to the next ch_attr
    if (! file.exists(predicted.values.dir)){next}
    
    ## Select the appropriate list of ecozones wrt ch_attr
    if (change.type == "always_treed") {
      ecozones <- params6$mapped.ecozones
    } else if (change.type == "fire") {
      ecozones <- params6$fire.ecozones
      boxplot.color <- params6$boxplots$color$fire
    } else if (change.type == "harvesting") {
      ecozones <- params6$harvest.ecozones
      boxplot.color <- params6$boxplots$color$harvesting
    }
    
    for (ez in 1:length(ecozones)) {
      
      ecozone.name <- names(ecozones[ez])
      ecozone.code <- ecozones[[ez]]
      
      print(sprintf("        Ecozone: %s", ecozone.name))
      
      ## Initialize empty list to store plots to then re-arrange them after: each ecozone with 2 boxplot time-series (+ count barplots)
      plots.list <- list()
      
      for (targ in params6$targ.names) {
          
        print(sprintf("            Targ: %s", targ))
        
        idx.targ <- which(targ == params6$targ.names)   ## numeric index useful to select elements of string vectors

        ## Set of target and ecozone-specific parameters
        if (ecozone.name %in% params3$sampled.ecozones.NONBOREAL) {
          ylim.targ.abs <- params6$ylim.targ$absolute[[targ]][["NONBOREAL"]]  ## larger range for y-axis in BC
        } else {
          ylim.targ.abs <- params6$ylim.targ$absolute[[targ]][["BOREAL"]]
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
                if (year != max(params6$mapped.years)) {
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
            fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("PixTimeSeries_1984_2016_%s_%s_%s_band5.pdf", temporal.data.type, change.type, gsub(" ", "", ecozone.name)), sep='')
            title.str <- sprintf("Band 5 %s time-series:\n%s, %s (%s pix.)", temporal.data.type, change.type, ecozone.name, params6$nr.pix.time.series)
            pdf(fig.name.str)
            theme_set(theme_gray(base_size = 18))
            time.series.band <- ggplot(temporal.traj.band.melted, aes(x=year, y=value, colour=variable, group=variable)) + 
              scale_colour_brewer(palette="Paired") +
              scale_x_continuous(breaks=seq(min(params6$mapped.years), max(params6$mapped.years), by=2), lim=c(min(params6$mapped.years), max(params6$mapped.years))) +
              xlab("year") +
              ylab("band 5") +
              coord_cartesian(ylim=params6$ylim.band5) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor.x=element_blank()) +
              geom_line() +
              ggtitle(title.str)
            print(time.series.band)
            dev.off()
            
          }  ## end if on targ == params6$targ.names[1] & method == params6$methods[1] to plot band 5
          
          
#### MELT FOREST ATTR DATA ---------------------------------------------------
          
          ## Build melted dt to plot time-series of boxplots
          
          ## Get landcover class for every pixel (rows) in each year (columns), used to determine trend type
          landcover.dt.subs <- data.table(pixID=yearly.predictions.full[, pixID])  ## initialize dt with pixID
          setkey(landcover.dt.subs, pixID)
          landcover.dt.subs <- merge(landcover.dt.subs, landcover, all=FALSE)  ## merge based on common pixID, as initial "metrics" files have been subset and contain thus less rows
          landcover.matrix <- data.matrix(landcover.dt.subs[, sprintf('%s', params6$mapped.years), with=FALSE])  ## corresponding matrix without pixID
          
          ## Matrices to get pre- and post-change majority label (coniferous/deciduous)
          landcover.class.pre.matrix <- landcover.matrix    ## initialize pre-distrubance land-cover matrix as a different matrix as landcover.class.pre.matrix will be subset when change.type != "always_treed"
          landcover.class.post.matrix <- landcover.matrix    ## same but for post-distrubance land-cover matrix

          if (change.type == "always_treed") {
            
            ## Years since change is a matrix of constant values (50), with some exceptions
            years.since.change.matrix <- matrix(rep(band.year.dt$YrsSince_GrCh, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE) 
            
            trend.type <- rep_len("alwTree", nrow(landcover.matrix))  ## if data is "always_treed", the only possible type of trend is "alwTree"
            
            xaxis.var <- "sampling_year"  ## dt column name for the variable to plot on the x-axis
            xlabel <- "year"  ## for no change x-axis represents the year 
  
          } else {
            
            ## Retrieve years since change data for every pixel (rows) in each year (columns)
            years.since.change.matrix <- matrix(nrow=nrow(band.year.dt), ncol=length(params6$mapped.years))  ## initialize empty matrix
            for (yr in 1:length(params6$mapped.years)) {
              year <- params6$mapped.years[yr]
              years.since.change.matrix[, yr] <- year - band.year.dt[, YrOfChange]
            }
            
            landcover.class.pre.matrix[years.since.change.matrix >= 0] <- 0   ## set to zero landcover values for years after change, including year of change
            landcover.class.post.matrix[years.since.change.matrix < 0] <- 0   ## set to zero landcover values for years before change, including year of change
            
            ## Determine pixel trend types
            trend.type <- apply(cbind(landcover.matrix, band.year.dt[, YrOfChange]), 1, determine_trend_types, treed.classes=params6$forest.classes$AllForest, no.tree.lag=params6$no.tree.lag, residual.lag=params6$residual.lag)
            ## TO BE USED TO CHECK
            ## determine_trend_types(cbind(landcover.matrix, band.year.dt[, YrOfChange])[14,], treed.classes=params6$forest.classes$AllForest, no.tree.lag=params6$no.tree.lag, residual.lag=params6$residual.lag)
            
            xaxis.var <- "year_since_ch"
            xlabel <- sprintf("years since %s", change.type)   ## for change x-axis represents the year since greatest change
            
          }
          
          ## Store counts per type of tend to se if we have too few samples in a plot
          list.trend.type.counts[[sprintf('%s_%s', ecozone.name, change.type)]]  <- table(trend.type)
          
          ## Get majority class in either pre- or post-change years
          forest.type.pre <- apply(landcover.class.pre.matrix, 1, maj_vote_non_zero)  # apply row-wise the majority vote function that returns most frequent label, excluding zeros 
          forest.type.pre.matrix <- matrix(rep(forest.type.pre, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE)  ## replicate vector columnwise to match shape of other matrices
          
          forest.type.post <- apply(landcover.class.post.matrix, 1, maj_vote_non_zero)   
          forest.type.post.matrix <- matrix(rep(forest.type.post, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE)
          
          ## Reshape trend types in a matrix
          trend.type.matrix <- matrix(rep(trend.type, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE)
          
          ## Stack together flattened matrices
          dt.yearly.pred.melted <- data.table(sampling_year=c( matrix(rep(params6$mapped.years, each=nrow(yearly.predictions.full)), nrow=nrow(yearly.predictions.full), byrow=FALSE) ),   ## c() reshapes matrix to vector
                                              pix_ID=c( matrix(rep(yearly.predictions.full$pixID, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE) ),
                                              year_since_ch=c( years.since.change.matrix ),    
                                              forest_attr=c( data.matrix(yearly.predictions.full[, sprintf('%s', params6$mapped.years), with=FALSE]) ),
                                              landcover_class=c( landcover.matrix ),  
                                              forest_type_pre=c( forest.type.pre.matrix ),
                                              forest_type_post=c( forest.type.post.matrix ),
                                              trend_type=c( trend.type.matrix )
                                  )
          
          ## Subset the dt to plot only attribute values of certain years (removal of "no tree lag" years). 
          ## Negative values of year_since_ch represent "time before change". 
          ## In the case of Always_treed this results in keeing all the rows, as time since change values are 50, which should always be > params6$start.change.plots 
          dt.yearly.pred.melted <- dt.yearly.pred.melted[(year_since_ch >= params6$start.change.plots & year_since_ch < 0) | 
                                                          year_since_ch >= params6$start.regrowth.ecoz.plots]  
          
          dt.yearly.pred.melted <- dt.yearly.pred.melted[landcover_class %in% params6$forest.classes$AllForest]  ## subset to plot only forested pixels
          dt.yearly.pred.melted[, year_since_ch:=as.factor(year_since_ch)]  ## set "year_since_ch" as factor to have one boxplot per year
          dt.yearly.pred.melted[, sampling_year:=as.factor(sampling_year)]  ## set "sampling_year" as factor to have one boxplot per year
          
#### BOXPLOTS ---------------------------------------------------
          
          ## Loop on forest classes (all forest, coniferous, deciduous)
          for (i in 1:length(params6$forest.classes)) {
              
            forest.class.name <- names(params6$forest.classes)[i]
            forest.labels <- params6$forest.classes[[i]]
            
            print(sprintf("                     Forest type: %s", forest.class.name))
            
            ## Subset to keep only pixels of a certain forest type based either on the pre- or post-change majority label
            if (params6$forest.type.period == "pre-change") {
              condition.forest.type <- "forest_type_pre %in% forest.labels" 
            } else if (params6$forest.type.period == "post-change") {
              condition.forest.type <- "forest_type_post %in% forest.labels"
            }
            
            for (t in 1:length(params6$trend.types)) {
              
              trend.type.name <- names(params6$trend.types)[t]
              trend.type.labels <- params6$trend.types[[t]]
            
              ## Only produce plots for stable treed pixels when trend.type.name == "AlwTreeInChange"
              if (change.type == "always_treed" & trend.type.name != "AlwTreeInChange") {
                next
              }
              
              ## Subset dt wrt forest and trend type
              condition.trend.type <- "trend_type %in% trend.type.labels"
              cmd <- sprintf("dt.yearly.pred.melted.to.plot <- dt.yearly.pred.melted[%s & %s]", condition.forest.type, condition.trend.type)
              eval(parse(text=cmd))
              
              ## Skip round of loop if the above conditions return an empty dt
              if (nrow(dt.yearly.pred.melted.to.plot) == 0) {
                next
              }
              
              ## Set different parameters/label types if we are plotting data for always_treed or regrow after distrubance (fire or harvesting)
              if (change.type == "always_treed") {
                trend.for.CAN.plot <- 'AlwTreeInChange'  ## type of data that has to participate to the Canada-wide plot (ecozones combined)
                plot.labels <- params6$mapped.years[seq(1, length(params6$mapped.years), by=2)]
                element_text.labels <- element_text(angle=45, hjust=1)
              } else {
                trend.for.CAN.plot <- 'Regrow'    ## type of data that has to participate to the Canada-wide plot (ecozones combined)
                plot.labels <- seq(params6$start.change.plots, max(params6$mapped.years)-min(params6$mapped.years)-1, by=3)
                element_text.labels <- element_text(angle=0)
              }
              
              ## When we have the trend type for the Canada-wide plot we grow a massive dt and we store it as a list by forest class, attribute (targ) and change type
              if (trend.type.name == trend.for.CAN.plot & temporal.data.type == 'fitted') {
                if (ez == 1) {   ## initialize dt for first ecozone...
                  list.dt.CAN.plot[[sprintf('%s_%s_%s', forest.class.name, targ, change.type)]] <- dt.yearly.pred.melted.to.plot[, c("sampling_year", "year_since_ch", "forest_attr"), with=FALSE]
                } else {  ## ...append for the rest
                  list.dt.CAN.plot[[sprintf('%s_%s_%s', forest.class.name, targ, change.type)]] <- rbind(list.dt.CAN.plot[[sprintf('%s_%s_%s',forest.class.name, targ, change.type)]], dt.yearly.pred.melted.to.plot[, c("sampling_year", "year_since_ch", "forest_attr"), with=FALSE])
                }
              }
              
              ## Store nr of pixels participating in each plot
              counts.per.year <- data.frame(table(dt.yearly.pred.melted.to.plot[, xaxis.var, with=FALSE]))  ## number of pixels in each yearly bin to plot reference histograms
              nr.plotted.pix <- length(unique(dt.yearly.pred.melted.to.plot$pix_ID))  ## number of pixels that participate to this plot
              list.nr.plotted.pix[[sprintf('%s_%s', forest.class.name, trend.type.name)]][ecozone.code, change.type] <- nr.plotted.pix

              ## Single code block to produce boxplots over time, either absolute time (year) or relative time (since change) for all change.type options (change and no-change classes)
              bplt <- ggplot(dt.yearly.pred.melted.to.plot) +
                    geom_boxplot(aes_string(x=xaxis.var, y="forest_attr"), outlier.shape=NA, notch=F) +
                    coord_cartesian(ylim=ylim.targ.abs) +
                    scale_x_discrete(drop=TRUE, breaks=plot.labels) +   ## to remove years from plot if they contain no data
                    ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
                    theme(axis.title.x=element_blank(), axis.text.x=element_text.labels, panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
              
              ## If "always_treed" run Mann-Kendall Test and add result as a plot subtitle
              if (change.type == "always_treed") {
                
                ## Mann-Kendall Test to see significance of trends in median values (Z > 0: increase in medians)
                medians <- ddply(dt.yearly.pred.melted.to.plot, .(sampling_year), summarize, med=median(forest_attr))
                if (!is.null(medians$med)) {
                  medians.ts <- ts(medians$med)
                  mk.test.res <- mk.test(medians.ts)
                  str <- sprintf("Mann-Kendall Test: Z = %2.2f, p-value = %2.3f", mk.test.res$statistic, mk.test.res$p.value)
                  bplt <- bplt + labs(subtitle=str)
                }
          
              ## otherwise if "Regrow" fit polynimials to quantiles in first 10 years
              } else if (trend.type.name == 'Regrow'){
              
                ## Compute quantiles by years_since_ch
                quantiles <- as.data.table( ddply(dt.yearly.pred.melted.to.plot, .(year_since_ch), summarize, 
                                                  c10=quantile(forest_attr, 0.1), 
                                                  c25=quantile(forest_attr, 0.25),
                                                  c50=quantile(forest_attr, 0.5), 
                                                  c75=quantile(forest_attr, 0.75),
                                                  c90=quantile(forest_attr, 0.9) ) )
                quantiles$year_since_ch <- as.numeric(as.character(quantiles$year_since_ch))
                
                ## Isolate quantile values to be fitted (zero plus the first 10 values after no tree lag)
                quantiles.to.fit <- quantiles[year_since_ch > params6$no.tree.lag & year_since_ch < params6$no.tree.lag+10]
                new.line <- data.frame(matrix(0, ncol=6, nrow=1))
                new.line[-1] <- rep(params6$start.val[[params6$targ.names.plots[idx.targ]]], 5)  ## add initial growth values for each quantile (dynamically changing based on targ name)
                names(new.line) <- names(quantiles.to.fit)
                quantiles.to.fit <- rbind(new.line, quantiles.to.fit)
                
                ## Define years since ch for which we want to predict a quantile value
                missing.data <- data.frame(year_since_ch=seq(0, params6$no.tree.lag))
                predicted.quantiles <- data.table(matrix(nrow = params6$no.tree.lag+1, ncol = 5))   ## initialize dt to store results
                names(predicted.quantiles) <- names(quantiles.to.fit)[-1]
                for (idx.quant in 1:ncol(predicted.quantiles)) {
                  quant.name <- names(predicted.quantiles)[idx.quant]
                  lm.model <- lm(data=quantiles.to.fit, sprintf('%s ~ poly(year_since_ch, 2)', quant.name))   ## fit a polynomial to each column of the dt with years_since_ch as x
                  # ggplot(quantiles.to.fit, aes(x=year_since_ch, y=c50)) +
                  #   geom_point() +
                  #   geom_smooth(method='lm', formula=y~poly(x, 2))
                  predicted.quantiles[[idx.quant]] <- predict(lm.model, newdata=missing.data)   ## predict a quantile value for each year we want a modeled value
                }
                
                ## To remove boxplot whiskers (weird extrapolated values) set same values for centiles 10 and 25 and for 90 and 75
                predicted.quantiles[, c10:=predicted.quantiles$c25]
                predicted.quantiles[, c90:=predicted.quantiles$c75]
                
                predicted.quantiles <- cbind(missing.data, predicted.quantiles)   ## add years_since_ch info
                
                ## Re-stack together the pre-change, the fitted and no-tree-lag quantile values
                boxplot.ts.data <- rbind(quantiles[year_since_ch<0], predicted.quantiles, quantiles[year_since_ch > params6$no.tree.lag])
                
                ## Add boolean indicator to plot extrapolated boxplots with different colors
                boxplot.ts.data[, extrapolation:=c(rep(F,nrow(quantiles[year_since_ch<0])), rep(T,nrow(predicted.quantiles)), rep(F, nrow(quantiles[year_since_ch > params6$no.tree.lag])))]
                
                boxplot.ts.data$year_since_ch <- as.factor(boxplot.ts.data$year_since_ch) ## to allow plotting with scale_x_discrete
                
                ## When trend type is 'Regrow' the boxplot with extrapolated values is overwritten on the same bplt ggplot object
                bplt <- ggplot(boxplot.ts.data, aes(x=year_since_ch, ymin=c10, lower=c25, middle=c50, upper=c75, ymax=c90, fill=extrapolation)) +
                  geom_boxplot(stat = "identity") +
                  scale_fill_manual(values = alpha(boxplot.color, c(1, params6$boxplots$alpha))) +
                  coord_cartesian(ylim=ylim.targ.abs) +
                  scale_x_discrete(drop=TRUE, breaks=plot.labels) +   ## to remove years from plot if they contain no data
                  ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
                  theme(axis.title.x=element_blank(), axis.text.x=element_text.labels, panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), legend.position="none")
              
              }  ## end else if (trend.type.name == 'Regrow')
                
              ## Store plot in dynamic list to be retrieved later on
              plots.list[[sprintf('Boxplot_%s_%s_%s_%s', params6$targ.names.plots[idx.targ], method, forest.class.name, trend.type.name)]] <- bplt
            
              ## Plot histogram of counts only for the first of the attributes, as it will be the same for the rest
              if (params6$targ.names.plots[idx.targ] == params6$targ.names.plots[1] & method == params6$methods[1]) {    
                
                ## Histogram with number of pixels in each time since change bin
                hst <- ggplot(counts.per.year, aes(x=Var1, y=Freq)) +
                      coord_fixed(ratio=1) +
                      geom_bar(stat="identity") +
                      ylab("pixel count") +
                      theme(aspect.ratio=0.25, axis.title.x=element_blank(), axis.text.x=element_text.labels, panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
                
                if (params6$boxplot.title == 'long') {
                  title.str <- sprintf("%s %s boxplots:\n%s, %s, %s, %s (%s pix.)", change.type, temporal.data.type, trend.type.name, forest.class.name, ecozone.name, method, nr.plotted.pix)
                } else if (params6$boxplot.title == 'short') {
                  title.str <- ecozone.name   ## if space between capital letters is needed sprintf("%s", gsub("([a-z])([A-Z])", "\\1 \\2", ecozone.name))
                }
                hst <- hst + ggtitle(title.str)
                
                plots.list[[sprintf('CountHist_%s_%s', forest.class.name, trend.type.name)]] <- hst
                
              }
            
            }   ## end for on params6$trend.types
            
          }  ## end for on params6$forest.classes
          
          
          #### To check values plotted #### 
          # quantile(yearly.predictions.full[, '2016', with=F][[1]], probs=c(0.1, 0.25, 0.5, 0.75, 0.9))  ## list subsetting of the one-column data.table gives a vector
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
            scale_x_continuous(breaks=seq(min(params6$mapped.years), max(params6$mapped.years), by=2), lim=c(min(params6$mapped.years), max(params6$mapped.years))) +
            xlab("year") +
            ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
            coord_cartesian(ylim=ylim.targ.abs) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
          
          
          ## Raw pixel time-series
          fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("PixTimeSeries_1984_2016_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
          title.str <- sprintf("%s %s time-series:\n%s, %s, %s (%s pix.)", params6$targ.names.plots[idx.targ], temporal.data.type, change.type, ecozone.name, method, params6$nr.pix.time.series)
          pdf(fig.name.str)
          theme_set(theme_gray(base_size = 18))
            time.series <- time.series.common + 
                            geom_line() + 
                            ggtitle(title.str)
          print(time.series)
          dev.off()
          
          ## Smoothed pixel time-series
          if (params6$plot.smoothed.ts) {
            fig.name.str <- file.path(FullPixelTimeSeries.subdir, sprintf("PixTimeSeries_1984_2016_smoothed_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, gsub(" ", "", ecozone.name), params6$targ.names.plots[idx.targ], method), sep='')
            title.str <- sprintf("%s %s smoothed time-series:\n%s, %s, %s (%s pix.)", params6$targ.names.plots[idx.targ], temporal.data.type, change.type, ecozone.name, method, params6$nr.pix.time.series)
            pdf(fig.name.str)
            theme_set(theme_gray(base_size = 18))
              time.series <- time.series.common + 
                              stat_smooth(method = "loess", formula = y ~ x, se = F) + 
                              ggtitle(title.str)
            print(time.series)
            dev.off()
          }
          
        } ## end for on params6$methods
      
      } ## end for on params6$targ.names.lg
    
      
      ## Retrieve plots from list of named plots and do figure with 3 sub-plots
      for (method in params6$methods) { 
        
        for (i in 1:length(params6$forest.classes)) {
          
          forest.class.name <- names(params6$forest.classes)[i]
          
          for (t in 1:length(params6$trend.types)) {
            
            trend.type.name <- names(params6$trend.types)[t]
            
            ## Only produce plots for stable treed pixels when trend.type.name == "AlwTreeInChange"
            if (change.type == "always_treed" & trend.type.name != "AlwTreeInChange") {
              next
            }
          
            ## Retrieve three sub-plots
            hst <- plots.list[[sprintf('CountHist_%s_%s', forest.class.name, trend.type.name)]]
            bplt.top <- plots.list[[sprintf('Boxplot_%s_%s_%s_%s', params6$targ.names.plots[1], method, forest.class.name, trend.type.name)]] 
            bplt.bottom <- plots.list[[sprintf('Boxplot_%s_%s_%s_%s', params6$targ.names.plots[2], method, forest.class.name, trend.type.name)]] 
            
            ## Skip round of loop if plots could not be produced due to empty dt 
            if (is.null(hst) | is.null(bplt.top) | is.null(bplt.bottom)) {
              next
            }
            
            ## Take plot title from histogram bc it is the valid for both cover_2m and elev_p95 subfigures
            bplt.top <- bplt.top + labs(title=hst$labels$title) 
            
            if (params6$boxplot.title == 'short') {
              bplt.top <- bplt.top + theme(plot.title=element_text(hjust=0.5, size=18, face="bold")) 
            }
            
            fig.name.str <- file.path(FullBoxPlots.subdir, sprintf("Boxplots_1984_2016_%s_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, forest.class.name, trend.type.name, gsub(" ", "", ecozone.name), method), sep='')
            pdf(fig.name.str)
            theme_set(theme_gray(base_size = 18))
              if (change.type == "always_treed") {
                  
                bplt.bottom <- bplt.bottom + theme(axis.title.x=element_text()) + xlab(xlabel)  ## receives here the xlabel, as either "years" or "years since <ch_attr>"
                
                bplt.top <- ggplotGrob(bplt.top)
                bplt.bottom <- ggplotGrob(bplt.bottom)
                grid.draw( rbind( bplt.top, bplt.bottom, size="last") )  ## to have plots neatly aligned under each other based on y-axis
                
              } else {
                
                bplt.top <- ggplotGrob(bplt.top)
                
                if (params6$plot.boxplot.counts) {
                  hst <- hst + theme(plot.title=element_blank(), axis.title.x=element_text()) + xlab(xlabel)  ## receives here the xlabel, as either "years" or "years since <ch_attr>"
                  hst <- ggplotGrob(hst)
                  bplt.bottom <- ggplotGrob(bplt.bottom)
                  grid.draw( rbind (bplt.top, bplt.bottom, hst, size="last") )
                } else {
                  bplt.bottom <- bplt.bottom + theme(plot.title=element_blank(), axis.title.x=element_text()) + xlab(xlabel)  ## receives here the xlabel, as either "years" or "years since <ch_attr>"
                  bplt.bottom <- ggplotGrob(bplt.bottom)
                  grid.draw( rbind (bplt.top, bplt.bottom, size="last") )
                }
                
              }
            dev.off()
          
          }
          
        }
        
      }
        
    } ## end for on params6$ecozone.code
    
  }  ## end for on params6a$change.types
  
}  ## end for on params6$temporal.data.types

## Plot fire vs harvest vs always_forest trends at CANADA level (all ecozones together)
for (i in 1:length(params6$forest.classes)) {
  
  forest.class.name <- names(params6$forest.classes)[i]
  
  line.plots <- list()
  for (targ in params6$targ.names) {
    
    idx.targ <- which(targ == params6$targ.names)   ## numeric index useful to select elements of string vectors
    
    ## Take respective Canada-wide dt
    dt.always_treed <- list.dt.CAN.plot[[sprintf('%s_%s_always_treed', forest.class.name, targ)]]
    dt.fire <- list.dt.CAN.plot[[sprintf('%s_%s_fire', forest.class.name, targ)]]
    dt.harvesting <- list.dt.CAN.plot[[sprintf('%s_%s_harvesting', forest.class.name, targ)]]
    
    ## Compute 3 quantiles only (Canada plots do not have 10th and 90th)
    quant.table.always_treed <- dt.always_treed[, (quantile(forest_attr, c(0.25, 0.5, 0.75)))]
    
    quant.table.fire <- dt.fire[, (quantile(forest_attr, c(0.25, 0.5, 0.75))), by=year_since_ch]
    quant.table.fire$quant <- rep(c("25thCentile", "50thCentile", "75thCentile"), length(unique(quant.table.fire$year_since_ch)))
    setnames(quant.table.fire, "V1", "fire")

    quant.table.harvesting <- dt.harvesting[, (quantile(forest_attr, c(0.25, 0.5, 0.75))), by=year_since_ch]
    quant.table.harvesting$quant <- rep(c("25thCentile", "50thCentile", "75thCentile"), length(unique(quant.table.harvesting$year_since_ch)))
    setnames(quant.table.harvesting, "V1", "harvesting")
    
    quant.table <- merge(quant.table.fire, quant.table.harvesting, by=c("year_since_ch", "quant"))
    quant.table[, always_treed:=rep(quant.table.always_treed, length(unique(quant.table$year_since_ch)))]
    
    quant.table[, year_since_ch:= as.numeric(as.character(quant.table$year_since_ch))]
    
    quant.table.wide <- cbind(data.table(year_since_ch=unique(quant.table$year_since_ch)),
                              quant.table[quant=="25thCentile", c("fire", "harvesting", "always_treed"), with=FALSE],
                              quant.table[quant=="50thCentile", c("fire", "harvesting", "always_treed"), with=FALSE],
                              quant.table[quant=="75thCentile", c("fire", "harvesting", "always_treed"), with=FALSE]
    )
    # quant.table.wide <- cbind(data.table(year_since_ch=unique(quant.table$year_since_ch)),
    #                           quant.table[quant=="25thCentile", c("fire", "harvesting"), with=FALSE],
    #                           quant.table[quant=="50thCentile", c("fire", "harvesting"), with=FALSE],
    #                           quant.table[quant=="75thCentile", c("fire", "harvesting"), with=FALSE]
    # )
    
    names(quant.table.wide) <- c("year_since_ch",
                                 "fire25th", "harvesting25th", "always_treed25th",
                                 "fire50th", "harvesting50th", "always_treed50th",
                                 "fire75th", "harvesting75th", "always_treed75th")
    # names(quant.table.wide) <- c("year_since_ch",
    #                              "fire25th", "harvesting25th",
    #                              "fire50th", "harvesting50th",
    #                              "fire75th", "harvesting75th")
    
    line.plots[[targ]] <- ggplot(quant.table.wide, aes(x=year_since_ch)) + 

                            geom_line(aes(y=fire25th), linetype="solid", size=0.5, color=params6$boxplots$color$fire) +
                            geom_line(aes(y=fire50th), linetype="solid", size=1, color=params6$boxplots$color$fire) +
                            geom_line(aes(y=fire75th), linetype="solid", size=0.5, color=params6$boxplots$color$fire) +
                            geom_ribbon(aes(ymin=fire25th, ymax=fire75th), alpha=0.1, fill=params6$boxplots$color$fire) +
      
                            geom_line(aes(y=harvesting25th), linetype="solid", size=0.5, color=params6$boxplots$color$harvesting) +
                            geom_line(aes(y=harvesting50th), linetype="solid", size=1, color=params6$boxplots$color$harvesting) +
                            geom_line(aes(y=harvesting75th), linetype="solid", size=0.5, color=params6$boxplots$color$harvesting) + 
                            geom_ribbon(aes(ymin=harvesting25th, ymax=harvesting75th), alpha=0.1, fill=params6$boxplots$color$harvesting) +
                            
                            geom_line(aes(y=always_treed25th), linetype="dashed", size=0.5, color=params6$boxplots$color$always.treed) +
                            geom_line(aes(y=always_treed50th), linetype="solid", size=1, color=params6$boxplots$color$always.treed) +
                            geom_line(aes(y=always_treed75th), linetype="dashed", size=0.5, color=params6$boxplots$color$always.treed) +
                            # geom_ribbon(aes(ymin=always_treed25th, ymax=always_treed75th), alpha=0.1, fill=params6$boxplots$color$always.treed) +

                            scale_x_continuous(limit=c(params6$start.regrowth.CAN.plots, max(quant.table$year_since_ch)), breaks=seq(params6$start.regrowth.CAN.plots, max(quant.table$year_since_ch), by=2)) +
                            ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
                            theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=0), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())

    
  }

  fig.name.str <- file.path(Temporal.subdir, "FullBoxPlots", sprintf("QuantilePlots_1984_2016_%s.pdf", forest.class.name), fsep = .Platform$file.sep)
  pdf(fig.name.str)
  theme_set(theme_gray(base_size = 18))
    bplt.bottom <- line.plots[[params6$targ.names[2]]] + theme(axis.title.x = element_text()) + xlab("years since disturbance")  ## receives here the xlabel, as either "years" or "years since <ch_attr>"
    bplt.top <- line.plots[[params6$targ.names[1]]]
    bplt.top <- ggplotGrob(bplt.top)
    bplt.bottom <- ggplotGrob(bplt.bottom)
    grid.draw( rbind( bplt.top, bplt.bottom, size="last") )  ## to have plots neatly aligned under each other based on y-axis
  dev.off()
  
}

list.trend.type.counts_path <- file.path(Temporal.subdir, 'trend_type_counts.Rdata', fsep = .Platform$file.sep) 
save(list.trend.type.counts, file=list.trend.type.counts_path)

list.nr.plotted.pix_path <- file.path(Temporal.subdir, 'nr_plotted_pix.Rdata', fsep = .Platform$file.sep) 
save(list.nr.plotted.pix, file=list.nr.plotted.pix_path)

## delete temp folder with temp rasters
unlink(raster_tmp_dir, recursive = T, force = T)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog6, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
