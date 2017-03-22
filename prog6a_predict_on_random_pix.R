## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog6a_predict_on_random_pix.R                           
## Objective: Predict forest attributes through time for randomly sampled pixels whose predictors are provided in a CSV file.

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - set appropriate FINAL_MODEL folder

# - Uncomment lines with ## ---------- UNCOMMENT ------------ 
# - ISSUE WITH BOREAL SHIELD WEST HAVING RANDOM TRENDS AND NOT STAIRS FOR YAI (see XXXXXXXXX)
# -

## SOLVED:
# -V yearly.predictors$YrsSince_GrCh is negative! -- added line to change values back to positive
# -V delete any pixID.file, pixID.by.eco.list reference -- removed this old trick to get same samples after subsetting both the datasets with predictors (with band 5 info) and predicted values: now we merge the datasets based on a common pixID "<UTMzone>_<index>"

#### READS/WRITES ------------------------------------------------------------

## READS:
# - 
## WRITES:
# - 

#### INIT --------------------------------------------------------------------

print('Prog6a: temporal predictions on random set of pixels') 

rm(list=ls())

param_file <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

param_file_prog3 = file.path(base_wkg_dir, 'AllUTMzones_params3.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog3)

stats_file_prog3 = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/FINAL_RUN/results/stats3.Rdata"
load(stats_file_prog3)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")

#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params6a <- list()

params6a$subsetting <- T  ## to subset the dataset to a nr of samples = params6a$nr.pts.subset (for debugging in development phase)
params6a$nr.pts.subset <- 2000
params6a$temporal.base.dir <- "E:/NTEMS/UTM_results/2017_02_27_SAMPLES_treed_changes_plots"  ## subdirectory to save files with predicted values
params6a$temporal.data.types <- c("fitted", "raw")   ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends
# params6a$temporal.data.types <- c("fitted")
params6a$change.types <-  c("always_treed", "fire", "harvesting")
# TODEL params6a$focus.on.UTM.zone <- F  ## to filter data to just a given UTM zone (test for fitted metrics to have a fair comparison with same pixels involved)
# TODEL params6a$UTM.zone <- '13S'   ## UTM zone to focus on
params6a$models.subdir <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/FINAL_RUN/Models"  ## model subdirectory to load models from
# TODEL params6a$subset.factor <- 4   ## to subset data to make it smaller, always applied bc we have too many samples
# params6a$mapped.years <- seq(from=1984, to=2012)
params6a$mapped.years <- seq(from=1984, to=1986)
# params6a$methods <- list("RF", "YAI")   ## has to be set as list("RF", "YAI"), that is both methods have to be run
params6a$methods <- list("YAI")   ## has to be set as list("RF", "YAI"), that is both methods have to be run
# params6a$mapped.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West",
#                                  "Hudson Plains", "Montane Cordillera", "Pacific Maritime", "Taiga Cordillera",
#                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")
# params6a$ecozone.codes <- c(2, 3, 4, 5, 6, 7, 9, 11, 15, 16, 17, 18)
# params6a$mapped.ecozones <- list("Boreal Plains", "Boreal Shield West", "Taiga Shield West")
# params6a$ecozone.codes <- c(4, 6, 18)
params6a$mapped.ecozones <- list("Boreal Shield West", "Boreal Plains")
params6a$ecozone.codes <- c(6, 4)
params6a$targ.names.temporal.lg <- c("elev_p95", "percentage_first_returns_above_2m")

param_file_prog6a = file.path(base_wkg_dir, 'AllUTMzones_params6a.Rdata', fsep = .Platform$file.sep)
save(params6a, file = param_file_prog6a)

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("latex2exp",    ## to have latex symbols in plots (not used in the end)
                      "xtable",   ## to print Latex tables
                      "rgeos",    ## spatial data handling...
                      "rgdal",
                      "sp",
                      "spdep",
                      "spatstat",
                      "ggplot2",
                      "gridExtra",
                      "GGally",
                      "psych",
                      "plyr",
                      "dplyr",   ## to be loaded before foreach to avoid "assertion failed" errors
                      "randomForest",
                      "rgl",
                      "yaImpute",
                      "caret",
                      "vegan",
                      "snow",
                      "lubridate", 
                      "doParallel", 
                      "foreach",
                      "data.table",
                      "RColorBrewer",
                      "gplots"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}


#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time


#### LOAD MODELS & READ DATA ------------------------------------------------

load(file.path(params6a$models.subdir, "YAI.Rdata", fsep=.Platform$file.sep))
load(file.path(params6a$models.subdir, "Ytrn.Rdata", fsep=.Platform$file.sep))

for (temporal.data.type in params6a$temporal.data.types) {
  
  print(sprintf("Temporal data: %s", temporal.data.type))

  for (change.type in params6a$change.types) {
    
    samples.dir <- file.path(params6a$temporal.base.dir, change.type, 'samples_u', fsep = .Platform$file.sep) 
    
    if (temporal.data.type == "fitted") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics_fitted', fsep = .Platform$file.sep)
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, 'predicted_values_fitted', fsep = .Platform$file.sep)
    } else if (temporal.data.type == "raw") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics', fsep = .Platform$file.sep) 
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, 'predicted_values', fsep = .Platform$file.sep) 
    }
    
    if (! file.exists(predicted.values.dir)){dir.create(predicted.values.dir, showWarnings = F, recursive = T)}
    
    nr.zero.samples.removed <- data.frame(matrix(nrow=length(params6a$mapped.ecozones), ncol=length(params6a$mapped.years)))
    rownames(nr.zero.samples.removed) <- params6a$mapped.ecozones
    colnames(nr.zero.samples.removed) <- paste0('yr', params6a$mapped.years)
    pct.aberr.elev <- nr.zero.samples.removed
    for (ez in 1:length(params6a$mapped.ecozones)) {
      
      ecozone.name <- params6a$mapped.ecozones[ez]
      ecozone.code <- params6a$ecozone.codes[ez]
      
      print(ecozone.name)
      
      file.name <- file.path(samples.dir, sprintf("%s_%s.csv", ecozone.code, ecozone.name), fsep = .Platform$file.sep)
      if (! file.exists(file.name)) {
        next
      }
      pix.IDs <- fread(file.name, header = F)
      
      predictions.list <- lapply(1:length(params6a$methods), function(x) NULL)
      names(predictions.list) <- params6a$methods
      
      for (yr in 1:length(params6a$mapped.years)) {
        
        if (yr %% 2 != 0) {
          print(params6a$mapped.years[yr])
        }
        
        ## Load dataset from CSV
        if (temporal.data.type == "fitted") {
          file.name <- file.path(predictors.dir, sprintf("%s_%s_%s_fitted.csv", ecozone.code, ecozone.name, params6a$mapped.years[yr]), fsep = .Platform$file.sep)
        } else {
          file.name <- file.path(predictors.dir, sprintf("%s_%s_%s.csv", ecozone.code, ecozone.name, params6a$mapped.years[yr]), fsep = .Platform$file.sep)
        }
        
        yearly.predictors <- fread(file.name, header = T)
        
        ## Change negative YrsSince_GrCh values back to positive
        yearly.predictors$YrsSince_GrCh <- abs(yearly.predictors$YrsSince_GrCh)
        
        pct.aberr.elev[ez, yr] <- sum(yearly.predictors$Elev > 5000) / nrow(yearly.predictors)
     
        ## Set unique pixel ID, with UTM zone info and index value
        yearly.predictors[, pixID:=sprintf('%s_%d', pix.IDs[[1]], pix.IDs[[2]])]
        yearly.predictors[, index:=NULL]

        
        # TODEL 
        ## Subset to UTM zone of focus (where we have both raw and fitted temporal data)
        # if (params6a$focus.on.UTM.zone){
        #   UTMz <- unlist(strsplit(yearly.predictors$pixID, '_'))[2*(1:length(yearly.predictors$pixID))-1]
        #   yearly.predictors <- yearly.predictors[UTMz==params6a$UTM.zone]
        # }
        # TODEL 
        
        # TODEL 
        # ## Subsample according to params6a$subset.factor and keep indices fixed per year
        # if (yr == 1) {
        #   set.seed(paramsGL$global.seed)
        #   idx.subsampling <- sample(c(T, F), size=nrow(yearly.predictors), prob=c(1/params6a$subset.factor, 1-(1/params6a$subset.factor)), replace=T)
        # }
        # yearly.predictors <- yearly.predictors[idx.subsampling, ]
        # TODEL 
        
        ## Preprocess columns
        yearly.predictors[, seq(1,ncol(yearly.predictors)-1):= na.roughfix(yearly.predictors[, !colnames(yearly.predictors) == 'pixID', with=FALSE])] ## substitute NAs with column median
        yearly.predictors[, Ch_attr:=factor(Ch_attr, levels=params3$Ch_attr.labels)]  ## set as factor with same levels as training set
        yearly.predictors[, Elev:=as.numeric(Elev)]

        samples.before <- nrow(yearly.predictors)
        yearly.predictors <- yearly.predictors[b1!=0 & b2!=0 & b3!=0 & b4!=0 & b5!=0 & b6!=0]
        nr.zero.samples.removed[ez, yr] <- samples.before - nrow(yearly.predictors)

        ## Trick to be used in the testing phase only to have smaller datasets and a faster debugging
        if (params6a$subsetting) {
          set.seed(paramsGL$global.seed)
          pix.idx <- sample(1:nrow(yearly.predictors), params6a$nr.pts.subset)
          yearly.predictors <- yearly.predictors[pix.idx, ]
        }

    #### RANDOM FOREST ---------------------------------------------------------------------

        if ("RF" %in% params6a$methods) {  ## run only if RF is in params6a$methods

          Y.map.predicted.RF <- data.table(rep(params6a$mapped.years[yr], nrow(yearly.predictors)))
          
          for (targ in params6a$targ.names.temporal.lg) {

              RF.model.path <- file.path(params6a$models.subdir, sprintf("RF_%s.Rdata", targ), fsep=.Platform$file.sep)
              load(RF.model.path)

              ## prediction on map pixels
              prediction.res <- predict(rf.RF, yearly.predictors[, stats3$final.predictors, with=FALSE], type="response", predict.all=T, nodes=F)

              Y.map.predicted.RF <- cbind(Y.map.predicted.RF, prediction.res$aggregate)

          }  ## end for on params6a$targ.names.temporal.lg

          setnames(Y.map.predicted.RF, c("Year", params6a$targ.names.temporal.lg))
          
          predictions.list$RF <- rbind(predictions.list$RF, list(Y.map.predicted.RF))
          

        }  ## end if on RF

    #### YAIMPUTE -------------------------------------------------------------------

        if ("YAI" %in% params6a$methods) {  ## run only if "YAI" is specified in methods

          Y.map.predicted.YAI <- data.table(rep(params6a$mapped.years[yr], nrow(yearly.predictors)))
          
          ids.val.predicted <- as.data.table( newtargets(yai.rf, yearly.predictors[, stats3$final.predictors, with=FALSE], k=1)["neiIdsTrgs"] )
          ids.val.predicted[, order:=seq(from=1, to=nrow(yearly.predictors[, stats3$final.predictors, with=FALSE]))]   ## add index telling the order of the samples before the merge
          colnames(ids.val.predicted)[1] <- "Predicted_FCID"
          Y.map.imputed <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
          setkey(Y.map.imputed, order)   ## after the merge resort with respect to that index...
          Y.map.imputed[, c("Predicted_FCID", "order"):=NULL]   ## ...and then remove it

          Y.map.predicted.YAI <- cbind(Y.map.predicted.YAI, Y.map.imputed)

          setnames(Y.map.predicted.YAI, c("Year", unlist(params3$targ.names.lg)))
          
          predictions.list$YAI <- rbind(predictions.list$YAI, list(Y.map.predicted.YAI))
          
        } ## end if on YAI
      

      }  ## end for on params6a$mapped.years
      
      ## Stitch together the list elements by column of lists
      melted.predictions <- lapply(predictions.list, function(x) rbindlist(x))
      
      ## Add pixel ID replicated by year to later reshape the data with one row per pixel
      melted.predictions <- lapply(melted.predictions, cbind, pixID=rep(yearly.predictors[,pixID], length(params6a$mapped.years)))
      
      lapply(seq_along(melted.predictions), function(idx) reshape_save(melted.predictions[[idx]], predicted.values.dir, change.type, names(melted.predictions)[idx], params6a$targ.names.temporal.lg, params6a$mapped.years, ecozone.code, ecozone.name))
      
    }  ## end for on params6a$mapped.ecozones
    
    ## ---------- UNCOMMENT ------------
    # zeroes.file <- file.path(params6a$temporal.base.dir, sprintf('ZeroSamplesRemovedByEcoByYear_%s.Rdata', temporal.data.type), fsep = .Platform$file.sep) 
    # save(nr.zero.samples.removed, file = zeroes.file)
    ## ---------- UNCOMMENT ------------
    
    aberr.elev.file <- file.path(params6a$temporal.base.dir, sprintf('AberrElevByEcoByYear_%s.csv', temporal.data.type), fsep = .Platform$file.sep) 
    write.csv(pct.aberr.elev, aberr.elev.file, row.names=T)
  
  }
  
}  ## end for on params6a$temporal.data.types


#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog6a, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)


