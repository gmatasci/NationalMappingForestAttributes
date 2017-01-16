## Project Name: NationalImputationForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog6a_predict_on_random_pix.R                           
## Objective: Predict forest attributes through time for randomly sampled pixels whose predictors are provided in a CSV file.

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# -

# - ISSUE WITH BOREAL SHIELD WEST HAVING RANDOM TRENDS AND NOT STAIRS FOR YAI

## SOLVED:
# -V yearly.predictors$YrsSince_GrCh is negative! -- added line to change values back to positive
# -V delete any pixID.file, pixID.by.eco.list reference -- removed this old trick to get same samples after subsetting both the datasets with predictors (with band 5 info) and predicted values: now we merge the datasets based on a common pixID "<UTMzone>_<index>"

#### READS/WRITES ------------------------------------------------------------

## READS:
# - 
## WRITES:
# - 

#### INIT --------------------------------------------------------------------

print('Prog4a: temporal predictions on random set of pixels') 

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

params4a <- list()

params4a$subsetting <- F  ## to subset the dataset to a nr of samples = params4a$nr.pts.subset (for debugging in development phase)
params4a$nr.pts.subset <- 2000
params4a$temporal.base.dir <- "E:/NTEMS/UTM_results/random_samples_always_forested"  ## subdirectory to save files with predicted values
params4a$temporal.data.types <- c("raw", "fitted")   ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends
# params4a$temporal.data.types <- c("fitted")
params4a$focus.on.UTM.zone <- T  ## to filter data to just a given UTM zone (test for fitted metrics to have a fair comparison with same pixels involved)
params4a$UTM.zone <- '13S'   ## UTM zone to focus on
params4a$models.subdir <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/FINAL_RUN/Models"  ## model subdirectory to load models from
params4a$subset.factor <- 4   ## to subset data to make it smaller, always applied bc we have too many samples
params4a$mapped.years <- seq(from=1984, to=2012)
# params4a$mapped.years <- seq(from=1984, to=1996)
params4a$methods <- list("RF", "YAI")   ## has to be set as list("RF", "YAI"), that is both methods have to be run
# params4a$mapped.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West",
#                                  "Hudson Plains", "Montane Cordillera", "Pacific Maritime", "Taiga Cordillera",
#                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")
# params4a$ecozone.codes <- c(2, 3, 4, 5, 6, 7, 9, 11, 15, 16, 17, 18)
params4a$mapped.ecozones <- list("Boreal Plains", "Boreal Shield West", "Taiga Shield West")
params4a$ecozone.codes <- c(4, 6, 18)
# params4a$mapped.ecozones <- list("Boreal Shield West", "Boreal Plains")
# params4a$ecozone.codes <- c(6, 4)
params4a$targ.names.temporal.lg <- c("elev_p95", "percentage_first_returns_above_2m")

param_file_prog4a = file.path(base_wkg_dir, 'AllUTMzones_params4a.Rdata', fsep = .Platform$file.sep)
save(params4a, file = param_file_prog4a)

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

load(file.path(params4a$models.subdir, "YAI.Rdata", fsep=.Platform$file.sep))
load(file.path(params4a$models.subdir, "Ytrn.Rdata", fsep=.Platform$file.sep))

for (temporal.data.type in params4a$temporal.data.types) {
  
  print(sprintf("Temporal data: %s", temporal.data.type))

  if (temporal.data.type == "fitted") {
    samples.dir <- file.path(params4a$temporal.base.dir, 'samples_fitted', fsep = .Platform$file.sep) 
    predictors.dir <- file.path(params4a$temporal.base.dir, 'metrics_fitted', fsep = .Platform$file.sep)
    predicted.values.dir <- file.path(params4a$temporal.base.dir, 'predicted_values_fitted', fsep = .Platform$file.sep)
  } else if (temporal.data.type == "raw") {
    samples.dir <- file.path(params4a$temporal.base.dir, 'samples', fsep = .Platform$file.sep) 
    predictors.dir <- file.path(params4a$temporal.base.dir, 'metrics', fsep = .Platform$file.sep) 
    predicted.values.dir <- file.path(params4a$temporal.base.dir, 'predicted_values', fsep = .Platform$file.sep) 
  }
  
  if (! file.exists(predicted.values.dir)){dir.create(predicted.values.dir, showWarnings = F, recursive = T)}
  
  nr.zero.samples.removed <- data.frame(matrix(nrow=length(params4a$mapped.ecozones), ncol=length(params4a$mapped.years)))
  rownames(nr.zero.samples.removed) <- params4a$mapped.ecozones
  colnames(nr.zero.samples.removed) <- paste0('yr', params4a$mapped.years)
  for (ez in 1:length(params4a$mapped.ecozones)) {
    
    ecozone.name <- params4a$mapped.ecozones[ez]
    ecozone.code <- params4a$ecozone.codes[ez]
    
    print(ecozone.name)
    
    file.name <- file.path(samples.dir, sprintf("%s_%s.csv", ecozone.code, ecozone.name), fsep = .Platform$file.sep)
    pix.IDs <- fread(file.name, header = F)
    
    melted.predictions <- NULL
  
    for (yr in 1:length(params4a$mapped.years)) {
      
      if (yr %% 2 != 0) {
        print(params4a$mapped.years[yr])
      }
      
      ## Load dataset from CSV
      file.name <- file.path(predictors.dir, sprintf("%s_%s_%s.csv", ecozone.code, ecozone.name, params4a$mapped.years[yr]), fsep = .Platform$file.sep)
      yearly.predictors <- fread(file.name, header = T)
      
      ## Change negative YrsSince_GrCh values back to positive
      yearly.predictors$YrsSince_GrCh <- abs(yearly.predictors$YrsSince_GrCh)
      
      ## Set unique pixel ID, with UTM zone info and index value
      yearly.predictors[, pixID:=sprintf('%s_%d', pix.IDs[[1]], pix.IDs[[2]])]
      yearly.predictors[, index:=NULL]
      
      ## Subset to UTM zone of focus (where we have both raw and fitted temporal data)
      if (params4a$focus.on.UTM.zone){
        UTMz <- unlist(strsplit(yearly.predictors$pixID, '_'))[2*(1:length(yearly.predictors$pixID))-1]   
        yearly.predictors <- yearly.predictors[UTMz==params4a$UTM.zone]
      }
      
      ## Subsample according to params4a$subset.factor and keep indices fixed per year
      if (yr == 1) {
        set.seed(paramsGL$global.seed)
        idx.subsampling <- sample(c(T, F), size=nrow(yearly.predictors), prob=c(1/params4a$subset.factor, 1-(1/params4a$subset.factor)), replace=T)
      }
      yearly.predictors <- yearly.predictors[idx.subsampling, ]
      
      ## Preprocess columns
      yearly.predictors[, seq(1,ncol(yearly.predictors)-1):= na.roughfix(yearly.predictors[, !colnames(yearly.predictors) == 'pixID', with=FALSE])] ## substitute NAs with column median
      yearly.predictors[, Ch_attr:=factor(Ch_attr, levels=params3$Ch_attr.labels)]  ## set as factor with same levels as training set
      yearly.predictors[, Elev:=as.numeric(Elev)]
      
      samples.before <- nrow(yearly.predictors)
      yearly.predictors <- yearly.predictors[b1!=0 & b2!=0 & b3!=0 & b4!=0 & b5!=0 & b6!=0]
      nr.zero.samples.removed[ez, yr] <- samples.before - nrow(yearly.predictors)
      
      ## Trick to be used in the testing phase only to have smaller datasets and a faster debugging
      if (params4a$subsetting) {
        set.seed(paramsGL$global.seed)
        pix.idx <- sample(1:nrow(yearly.predictors), params4a$nr.pts.subset)
        yearly.predictors <- yearly.predictors[pix.idx, ]
      }
      
      Y.map.predicted.RF <- Y.map.predicted.YAI <- data.table(rep(params4a$mapped.years[yr], nrow(yearly.predictors)))
      
  #### RANDOM FOREST ---------------------------------------------------------------------
      
      if ("RF" %in% params4a$methods) {  ## run only if RF is in params4a$methods
  
        for (targ in params4a$targ.names.temporal.lg) {
  
            RF.model.path <- file.path(params4a$models.subdir, sprintf("RF_%s.Rdata", targ), fsep=.Platform$file.sep)
            load(RF.model.path)
  
            ## prediction on map pixels
            prediction.res <- predict(rf.RF, yearly.predictors[, stats3$final.predictors, with=FALSE], type="response", predict.all=T, nodes=F)
  
            Y.map.predicted.RF <- cbind(Y.map.predicted.RF, prediction.res$aggregate)
  
        }  ## end for on params4a$targ.names.temporal.lg
  
        setnames(Y.map.predicted.RF, c("Year", params4a$targ.names.temporal.lg))
  
      }  ## end if on RF
          
  #### YAIMPUTE -------------------------------------------------------------------
          
      if ("YAI" %in% params4a$methods) {  ## run only if "YAI" is specified in methods
  
        ids.val.predicted <- as.data.table( newtargets(yai.rf, yearly.predictors[, stats3$final.predictors, with=FALSE], k=1)["neiIdsTrgs"] )
        ids.val.predicted[, idx:=seq(from=1, to=nrow(yearly.predictors[, stats3$final.predictors, with=FALSE]))]   ## add index telling the order of the samples before the merge
        colnames(ids.val.predicted)[1] <- "Predicted_FCID"
        Y.map.imputed <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
        setkey(Y.map.imputed, idx)   ## after the merge resort with respect to that index...
        Y.map.imputed[, c("Predicted_FCID", "idx"):=NULL]   ## ...and then remove it
  
        Y.map.predicted.YAI <- cbind(Y.map.predicted.YAI, Y.map.imputed)
  
      } ## end if on YAI
  
      setnames(Y.map.predicted.YAI, c("Year", unlist(params3$targ.names.lg)))
  
      ## Grow list of lists (nr. years x 2)
      melted.predictions <- rbind(melted.predictions, list(Y.map.predicted.RF, Y.map.predicted.YAI))
      
      
    }  ## end for on params4a$mapped.years
    
    ## Stitch together the list elements by column of lists
    melted.predictions.RF <- rbindlist(melted.predictions[,1])  ## first is for RF
    melted.predictions.YAI <- rbindlist(melted.predictions[,2])   ## second is for YAI
  
    ## Add pixel ID replicated by year to later reshape the data with one row per pixel
    melted.predictions.RF[, pixID:=rep(yearly.predictors[,pixID], length(params4a$mapped.years))]
    melted.predictions.YAI[, pixID:=rep(yearly.predictors[,pixID], length(params4a$mapped.years))]
  
    ## Loop to select only the targets of interest
    for (targ in params4a$targ.names.temporal.lg) {
  
      
      # XXXXXXXXXXXXX
      # XXXXXXXXXXXXX
      # Check if all Ok with stairs in predictions.RF and predictions.YAI for both ecozones, now that "Boreal Shield West" is the 1st one, and check if errors are now for 2nd ecozone

      
      
      ## Reshape RF data to have one row per pixel (pixID) and one column per year and then add ecozone ID
      predictions.RF <- reshape(melted.predictions.RF[, c("pixID", "Year", targ), with=FALSE], idvar = "pixID", timevar = "Year", direction = "wide")
      setnames(predictions.RF, c('pixID', as.character(params4a$mapped.years)))
      predictions.RF[, ecozone:=ecozone.code]
  
      ## Save in a file with name varying wrt target and ecozone
      file.name.melted.predictions.RF <- file.path(predicted.values.dir, sprintf("persistant_treed_RF_%s_%s.csv", targ, gsub(' ', '', ecozone.name)), fsep = .Platform$file.sep)
      fwrite(predictions.RF, file.name.melted.predictions.RF)
  
      ## Same for YAI
      predictions.YAI <- reshape(melted.predictions.YAI[, c("pixID", "Year", targ), with=FALSE], idvar = "pixID", timevar = "Year", direction = "wide")
      setnames(predictions.YAI, c('pixID', as.character(params4a$mapped.years)))
      predictions.YAI[, ecozone:=ecozone.code]
  
      file.name.melted.predictions.YAI <- file.path(predicted.values.dir, sprintf("persistant_treed_YAI_%s_%s.csv", targ, gsub(' ', '', ecozone.name)), fsep = .Platform$file.sep)
      fwrite(predictions.YAI, file.name.melted.predictions.YAI)
  
    }
    
  }  ## end for on params4a$mapped.ecozones
  
  zeroes.file <- file.path(params4a$temporal.base.dir, sprintf('ZeroSamplesRemovedByEcoByYear_%s.Rdata', temporal.data.type), fsep = .Platform$file.sep) 
  save(nr.zero.samples.removed, file = zeroes.file)

}  ## end for on params4a$temporal.data.types


#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog4a, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)


