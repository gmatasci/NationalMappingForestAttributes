## Project Name: NationalImputationForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog6b_predict_on_val_pix.R                          
## Objective: Predict forest attributes through time for validation pixels whose predictors are provided in a CSV file.

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:

# - add prediction/imput on csv with val points 

## SOLVED:
# -

#### READS/WRITES ------------------------------------------------------------

## READS:
# - 
## WRITES:
# - 

#### INIT --------------------------------------------------------------------

start.message <- sprintf("Prog4b: temporal predictions on validation set, started running on %s", Sys.time())
print(start.message)

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

## Update results directory to point to the FINAL results one
base_results_dir <- file.path(base_dir, "FINAL_RUN", "results", fsep=.Platform$file.sep)

## subdirectory to save files with predicted values
temporal.base.dir <- "E:/NTEMS/UTM_results/samples_plots"
predictors.dir <- file.path(temporal.base.dir, 'metrics', fsep = .Platform$file.sep) 
predicted.values.dir <- file.path(temporal.base.dir, 'predicted_values', fsep = .Platform$file.sep) 
if (! file.exists(predicted.values.dir)){dir.create(predicted.values.dir, showWarnings = F, recursive = T)}

## model subdirectory to load models from
models.subdir <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/FINAL_RUN/Models"


#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params4b <- list()

## Actual parameters to be used
params4b$subsetting <- F  ## to subset the dataset to a nr of samples = params3$nr.pts.plot (for debugging in development phase)
params4b$nr.pts.subset <- 2000
params4b$mapped.years <- seq(from=1984, to=2012)
# params4b$mapped.years <- seq(from=1984, to=1986)
params4b$methods <- c("RF", "YAI") ## accepts "RF" and/or "YAI", used to run analyses only for the specified methods

params4b$targ.names.temporal.lg <- c("elev_p95", "percentage_first_returns_above_2m")

param_file_prog4b = file.path(base_wkg_dir, 'AllUTMzones_params4b.Rdata', fsep = .Platform$file.sep) 
save(params4b, file = param_file_prog4b)

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

#### LOAD MODELS & READ DATA --------------------------------------------------------------

load(file.path(models.subdir, "YAI.Rdata", fsep=.Platform$file.sep))
load(file.path(models.subdir, "Ytrn.Rdata", fsep=.Platform$file.sep))
  
plots.info <- fread("D:/Research/ANALYSES/NationalMappingForestAttributes/Rcode_supercomputer_Txomin/Get_plots_coords/plot_UTM_info.csv", header = T)

melted.predictions <- NULL

for (yr in 1:length(params4b$mapped.years)) {
  
  print(params4b$mapped.years[yr])
  
  ## Load dataset from CSV
  file.name <- file.path(predictors.dir, sprintf("0_plot_UTM_indices_%s.csv", params4b$mapped.years[yr]), fsep = .Platform$file.sep)
  yearly.predictors <- fread(file.name, header = T)
  
  ## Keep only validation data
  yearly.predictors <- yearly.predictors[plots.info$TV=="VALIDATION", ] 
  
  ## Preprocess rows and columns
  yearly.predictors <- na.roughfix(yearly.predictors) ## substitute NAs with column median
  yearly.predictors[, Ch_attr:=factor(Ch_attr, levels=params3$Ch_attr.labels)]  ## set as factor with same levels as training set
  yearly.predictors[, Elev:=as.numeric(Elev)]  
  yearly.predictors[, YrsSince_GrCh:=as.numeric(YrsSince_GrCh)]  
  
  ## Trick to be used in the testing phase only to have smaller datasets and a faster debugging
  if (params4b$subsetting) {
    set.seed(paramsGL$global.seed)
    pix.idx <- sample(1:nrow(yearly.predictors), params4b$nr.pts.subset)
    yearly.predictors <- yearly.predictors[pix.idx, ]
  }
  
  Y.map.predicted.RF <- Y.map.predicted.YAI <- data.table(rep(params4b$mapped.years[yr], nrow(yearly.predictors)))
  
  
#### RANDOM FOREST ---------------------------------------------------------------------
  
  if ("RF" %in% params4b$methods) {  ## run only if RF is in params4b$methods
    
    for (targ in params4b$targ.names.temporal.lg) {
    
        RF.model.path <- file.path(models.subdir, sprintf("RF_%s.Rdata", targ), fsep=.Platform$file.sep)
        load(RF.model.path)
        
        ## prediction on map pixels
        prediction.res <- predict(rf.RF, yearly.predictors[, stats3$final.predictors, with=FALSE], type="response", predict.all=T, nodes=F)
        
        Y.map.predicted.RF <- cbind(Y.map.predicted.RF, prediction.res$aggregate)
        
    }  ## end for on params4b$targ.names.temporal.lg
    
    setnames(Y.map.predicted.RF, c("Year", params4b$targ.names.temporal.lg))
    
  }  ## end if on RF
      
#### YAIMPUTE -------------------------------------------------------------------
      
  if ("YAI" %in% params4b$methods) {  ## run only if "YAI" is specified in methods
    
    ids.val.predicted <- as.data.table( newtargets(yai.rf, yearly.predictors[, stats3$final.predictors, with=FALSE], k=1)["neiIdsTrgs"] )
    ids.val.predicted[, idx:=seq(from=1, to=nrow(yearly.predictors[, stats3$final.predictors, with=FALSE]))]
    colnames(ids.val.predicted)[1] <- "Predicted_FCID"
    Y.map.imputed <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
    setkey(Y.map.imputed, idx)
    Y.map.imputed[, c("Predicted_FCID", "idx"):=NULL]
    
    Y.map.predicted.YAI <- cbind(Y.map.predicted.YAI, Y.map.imputed)
    
  } ## end if on YAI

  setnames(Y.map.predicted.YAI, c("Year", unlist(params3$targ.names.lg)))
  
  ## Grow list of lists (nr. years x 2)
  melted.predictions <- rbind(melted.predictions, list(Y.map.predicted.RF, Y.map.predicted.YAI))
  
  
}  ## end for on params4b$mapped.years

## Stitch together the list elements by column of lists
melted.predictions.RF <- rbindlist(melted.predictions[,1])  ## first is for RF
melted.predictions.YAI <- rbindlist(melted.predictions[,2])   ## second is for YAI

## Add pixel ID replicated by year to later reshape the data with one row per pixel
melted.predictions.RF[, pixID:=rep(1:sum(melted.predictions.RF$Year==params4b$mapped.years[1]), length(params4b$mapped.years))]
melted.predictions.YAI[, pixID:=rep(1:sum(melted.predictions.YAI$Year==params4b$mapped.years[1]), length(params4b$mapped.years))]

## Loop to select only the targets of interest
for (targ in params4b$targ.names.temporal.lg) {
  
  ## Reshape RF data to have one row per pixel (pixID) and one column per year and then add ecozone ID
  predictions.RF <- reshape(melted.predictions.RF[, c("pixID", "Year", targ), with=FALSE], idvar = "pixID", timevar = "Year", direction = "wide")
  predictions.RF[, pixID:=NULL]
  setnames(predictions.RF, as.character(params4b$mapped.years))

  ## Save in a file with name varying wrt target and ecozone
  file.name.melted.predictions.RF <- file.path(predicted.values.dir, sprintf("time_series_val_set_RF_%s.csv", targ), fsep = .Platform$file.sep)
  fwrite(predictions.RF, file.name.melted.predictions.RF)
  
  ## Same for YAI
  predictions.YAI <- reshape(melted.predictions.YAI[, c("pixID", "Year", targ), with=FALSE], idvar = "pixID", timevar = "Year", direction = "wide")
  predictions.YAI[, pixID:=NULL]
  setnames(predictions.YAI, as.character(params4b$mapped.years))

  file.name.melted.predictions.YAI <- file.path(predicted.values.dir, sprintf("time_series_val_set_YAI_%s.csv", targ), fsep = .Platform$file.sep)
  fwrite(predictions.YAI, file.name.melted.predictions.YAI)
  
}
  
#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog4b, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)

