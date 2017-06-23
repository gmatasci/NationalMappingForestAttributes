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

# params6$experiment.name <- "fitted"  ## matching to folder names in E:\NTEMS\UTM_results\2017_02_27_SAMPLES_treed_changes_plots\<ch_attr> and in D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes
params6$experiment.name <- "fitted_NO_CHATTR"
# params6$experiment.name <- "fitted_NO_TEMP"
# params6$experiment.name <- "fitted_NO_CHATTR_LAT_LONG"

# params6$temporal.data.types <- c("raw", "fitted")  ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends
params6$temporal.data.types <- c("fitted")

params6$change.types <-  c("always_treed", "fire", "harvesting")
# params6$change.types <- c("always_treed", "harvesting")
# params6$change.types <-  c("harvesting")

params6$methods <- c("YAI")   ## for which methods to run the analysis, can accept both c("RF", "YAI") 

params6$mapped.years <- seq(from=1984, to=2012)     ## only used for val boxplots when reading back the images
# params6$mapped.years <- seq(from=1984, to=1988)
# params6$mapped.years <- seq(from=2008, to=2012)

# params6$forest.type.period <- "pre-change"   ## define forest type based on either pre- or post-change majority class
params6$forest.type.period <- "post-change"  

params6$forest.classes <- list("AllForest"=c(81, 210, 220, 230), "WetlandTreed"=81, "ConiferousForest"=210, "DeciduousForest"=220, "MixedForest"=230)

params6$boxplot.title <- 'short'  ## set to 'long' for testing phase to have a telling text in the figure, set to 'short' when producing figures for paper

# params6$start.change.plots <- -50  ## year at which to start plotting the boxplot time-series (x-axis lower limit)
params6$start.change.plots <- -3

params6$plot.smoothed.ts <- F    ## whether to produce Loess smoothed plots of single pixel time-series

params6$mapped.ecozones <- list("Atlantic Maritime"=2, "Boreal Cordillera"=3, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Hudson Plains"=7, "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Cordillera"=15,
                                 "Taiga Plains"=16, "Taiga Shield East"=17, "Taiga Shield West"=18)

params6$harvest.ecozones <- list("Atlantic Maritime"=2, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Plains"=16)

params6$fire.ecozones <- list("Boreal Cordillera"=3, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Hudson Plains"=7, "Montane Cordillera"=9, "Taiga Cordillera"=15,
                                 "Taiga Plains"=16, "Taiga Shield East"=17, "Taiga Shield West"=18)

# params6$mapped.ecozones <- list("Pacific Maritime"=11)

# params6$harvest.ecozones <- list("Atlantic Maritime"=2)

# params6$fire.ecozones <- list("Boreal Cordillera"=3)


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

for (temporal.data.type in params6$temporal.data.types) {
  
  print(sprintf("Temporal data: %s", temporal.data.type))
  
  for (change.type in params6$change.types) {
    
    print(sprintf("    Ch_attr: %s", change.type))
    
    samples.dir <- file.path(params6a$temporal.base.dir, change.type, 'samples_u', fsep = .Platform$file.sep)
    landcover.dir <- file.path(params6a$temporal.base.dir, change.type, 'landcover', fsep = .Platform$file.sep)
    
    if (temporal.data.type == "fitted") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics_fitted', fsep = .Platform$file.sep)
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, sprintf('predicted_values_%s', params6$experiment.name), fsep = .Platform$file.sep)
      Temporal.subdir <- file.path(sprintf('%s_%s', base_figures_dir, params6$experiment.name), "TemporalAnalyses", fsep = .Platform$file.sep)
    } else if (temporal.data.type == "raw") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics', fsep = .Platform$file.sep) 
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, 'predicted_values', fsep = .Platform$file.sep) 
      Temporal.subdir <- file.path(base_figures_dir, "TemporalAnalyses", fsep = .Platform$file.sep)
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
    } else if (change.type == "harvesting") {
      ecozones <- params6$harvest.ecozones
    }
    
    for (ez in 1:length(ecozones)) {
      
      ecozone.name <- names(ecozones[ez])
      ecozone.code <- ecozones[[ez]]
      
      print(sprintf("        Ecozone: %s", ecozone.name))
      
      plots.list <- list()
      
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
              theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor.x=element_blank()) +
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
          landcover.class.pre.matrix <- landcover.matrix    ## store as a different matrix as landcover.class.pre.matrix will be subset when change.type != "always_treed"
          landcover.class.post.matrix <- landcover.matrix    ## store as a different matrix as landcover.class.post.matrix will be subset when change.type != "always_treed"
          
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
            
            landcover.class.pre.matrix[years.since.change.matrix >= 0] <- 0   ## set to zero landcover values for years after change, including year of change
            landcover.class.post.matrix[years.since.change.matrix < 0] <- 0   ## set to zero landcover values for years after change, including year of change
            
            xaxis.var <- "year_since_ch"
            xlabel <- sprintf("years since %s", change.type)   ## for change x-axis represents the year since greatest change
            
          }
          
          ## Get majority class in either pre- or post-change years
          forest.type.pre <- apply(landcover.class.pre.matrix, 1, maj_vote_non_zero)  # apply row-wise the majority vote function that returns most frequent label, excluding zeros 
          forest.type.pre.matrix <- matrix(rep(forest.type.pre, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE)  ## replicate vector columnwise to match shape of other matrices
          
          forest.type.post <- apply(landcover.class.post.matrix, 1, maj_vote_non_zero)   
          forest.type.post.matrix <- matrix(rep(forest.type.post, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE)
          
          
          ## Stack together flattened matrices
          dt.yearly.pred.melted <- data.table(sampling_year=c( matrix(rep(params6$mapped.years, each=nrow(yearly.predictions.full)), nrow=nrow(yearly.predictions.full), byrow=FALSE) ),   ## c() reshapes matrix to vector
                                              pix_ID=c( matrix(rep(yearly.predictions.full$pixID, each=length(params6$mapped.years)), ncol=length(params6$mapped.years), byrow=TRUE) ),
                                              year_since_ch=c( years.since.change.matrix ),    
                                              forest_attr=c( data.matrix(yearly.predictions.full[, sprintf('%s', params6$mapped.years), with=FALSE]) ),
                                              landcover_class=c( landcover.matrix ),  
                                              forest_type_pre=c( forest.type.pre.matrix ),
                                              forest_type_post=c( forest.type.post.matrix )
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
            
            ## Subset to keep only pixels of a certain forest type based either on the pre- or post-change majority label
            if (params6$forest.type.period == "pre-change") {
              dt.yearly.pred.melted.to.plot <- dt.yearly.pred.melted[forest_type_pre %in% forest.labels] 
            } else if (params6$forest.type.period == "post-change") {
              dt.yearly.pred.melted.to.plot <- dt.yearly.pred.melted[forest_type_post %in% forest.labels]  
            }
            
            counts.per.year <- data.frame(table(dt.yearly.pred.melted.to.plot[, xaxis.var, with=FALSE]))  ## number of pixels in each yearly bin to plot reference histograms
            nr.plotted.pix <- length(unique(dt.yearly.pred.melted.to.plot$pix_ID))  ## number of pixels that participate to this plot
            
            ## Single code block to produce boxplots over time, either absolute time (year) or relative time (since change) for all change.type options (change and no-change classes)
            bplt <- ggplot(dt.yearly.pred.melted.to.plot) +
                  geom_boxplot(aes_string(x=xaxis.var, y="forest_attr"), outlier.shape=NA, notch=F) +
                  coord_cartesian(ylim=ylim.targ.abs) +
                  ylab(sprintf("%s [%s]", params6$targ.names.plots[idx.targ], params6$targ.units[idx.targ])) +
                  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
            
            if (change.type == "always_treed") {
              
              ## Mann-Kendall Test to see significance of trends in median values (Z > 0: increase in medians)
              medians <- ddply(dt.yearly.pred.melted.to.plot, .(sampling_year), summarize, med=median(forest_attr))
              if (!is.null(medians$med)) {
                medians.ts <- ts(medians$med)
                mk.test.res <- mk.test(medians.ts)
                str <- sprintf("Mann-Kendall Test: Z = %2.2f, p-value = %2.3f", mk.test.res$Z, mk.test.res$pvalue)
                bplt <- bplt + labs(subtitle=str)
              }
        
            }
            
            plots.list[[sprintf('Boxplot_%s_%s_%s', params6$targ.names.plots[idx.targ], method, forest.class.name)]] <- bplt  ## save plot in list to be retrieved later on
          
            if (params6$targ.names.plots[idx.targ] == params6$targ.names.plots[1]) {    ## save bin counts only for the first of the attributes, as it will be the same for the rest
            
              ## Histogram with number of pixels in each time since change bin
              hst <- ggplot(counts.per.year, aes(x=Var1, y=Freq)) +
                    coord_fixed(ratio=1) +
                    geom_bar(stat="identity") +
                    ylab("pixel count") +
                    theme(aspect.ratio=0.25, axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
              
              if (params6$boxplot.title == 'long') {
                title.str <- sprintf("%s boxplots:\n%s, %s, %s, %s (%s pix.)", temporal.data.type, change.type, forest.class.name, ecozone.name, method, nr.plotted.pix)
              } else if (params6$boxplot.title == 'short') {
                title.str <- ecozone.name   ## if space between capital letters is needed sprintf("%s", gsub("([a-z])([A-Z])", "\\1 \\2", ecozone.name))
              }
              hst <- hst + ggtitle(title.str)
              
              plots.list[[sprintf('CountHist_%s', forest.class.name)]] <- hst
              
            }
            
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
            theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
          
          
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
    
      
      ## Retrieve plots from list of named plots and do figure with 3 sub-plots
      for (method in params6$methods) { 
        
        for (i in 1:length(params6$forest.classes)) {
          
          forest.class.name <- names(params6$forest.classes)[i]
          
          hst <- plots.list[[sprintf('CountHist_%s', forest.class.name)]]
          bplt.elev_p95 <- plots.list[[sprintf('Boxplot_elev_p95_%s_%s', method, forest.class.name)]] 
          bplt.cover_2m <- plots.list[[sprintf('Boxplot_cover_2m_%s_%s', method, forest.class.name)]] 
          
          bplt.elev_p95 <- bplt.elev_p95 + labs(title=hst$labels$title) + theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
          
          fig.name.str <- file.path(FullBoxPlots.subdir, sprintf("Boxplots_1984_2012_%s_%s_%s_%s_%s.pdf", temporal.data.type, change.type, forest.class.name, gsub(" ", "", ecozone.name), method), sep='')
          pdf(fig.name.str)
          theme_set(theme_gray(base_size = 15))
            if (change.type == "always_treed") {
                
              bplt.cover_2m <- bplt.cover_2m + theme(axis.title.x = element_text(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab(xlabel)  ## receives here the xlabel, as either "years" or "years since <ch_attr>"
              
              bplt.elev_p95 <- ggplotGrob(bplt.elev_p95)
              bplt.cover_2m <- ggplotGrob(bplt.cover_2m)
              grid.draw( rbind( bplt.elev_p95, bplt.cover_2m ) )  ## to have plots neatly aligned under each other based on y-axis
              
            } else {
              
              hst <- hst + theme(plot.title = element_blank(), axis.title.x = element_text(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab(xlabel)  ## receives here the xlabel, as either "years" or "years since <ch_attr>"
              hst <- ggplotGrob(hst)
              bplt.elev_p95 <- ggplotGrob(bplt.elev_p95)
              bplt.cover_2m <- ggplotGrob(bplt.cover_2m)
              grid.draw( rbind (bplt.elev_p95, bplt.cover_2m, hst) )

            }
          dev.off()
          
        }
        
      }
        
    } ## end for on params6$ecozone.code
    
  }  ## end for on params6a$change.types
  
}  ## end for on params6$temporal.data.types


## delete temp folder with temp rasters
unlink(raster_tmp_dir, recursive = T, force = T)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog6, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
