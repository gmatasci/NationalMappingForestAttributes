#### CODE INFOS -------------------------------------------------------------------

## Project Name: NationalImputationForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name:                            
## Objective: Descriptive statistics and plots, variable selection and Random Forest prediction/imputation

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:



## SOLVED:
# -

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

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatImp.R")

## subdirectory to save files with predicted values
temporal.base.dir <- "E:/NTEMS/UTM_results/random_samples_always_forested"
predictors.dir <- file.path(temporal.base.dir, 'metrics', fsep = .Platform$file.sep) 
predicted.values.dir <- file.path(temporal.base.dir, 'predicted_values', fsep = .Platform$file.sep) 
if (! file.exists(predicted.values.dir)){dir.create(predicted.values.dir, showWarnings = F, recursive = T)}

## model subdirectory to load models from
models.subdir <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/FINAL_RUN/Models"


#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params4a <- list()

## Actual parameters to be used
params4a$subsetting <- T  ## to subset the dataset to a nr of samples = params3$nr.pts.plot (for debugging in development phase)
# params4a$mapped.years <- seq(from=1984, to=2012)
params4a$mapped.years <- seq(from=1984, to=1986)
params4a$methods <- list("RF", "YAI") ## accepts "RF" and/or "YAI", used to run analyses only for the specified methods
# params4a$methods <- list("YAI") ## accepts "RF" and/or "YAI", used to run analyses only for the specified methods

# params4a$mapped.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", 
#                                  "Hudson Plains", "Montane Cordillera", "Pacific Maritime", "Taiga Cordillera", 
#                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")
# params4a$ecozone.codes <- c(2, 3, 4, 5, 6, 7, 9, 11, 15, 16, 17, 18)

params4a$mapped.ecozones <- c("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains")
params4a$ecozone.codes <- c(2, 3, 4)

params4a$nr.pts.subset <- 2000

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

#### LOAD MODELS & READ DATA --------------------------------------------------------------

load(file.path(models.subdir, "YAI.Rdata", fsep=.Platform$file.sep))
load(file.path(models.subdir, "Ytrn.Rdata", fsep=.Platform$file.sep))

for (ez in 1:length(params4a$mapped.ecozones)) {
  
  ecozone.name <- params4a$mapped.ecozones[ez]
  ecozone.code <- params4a$ecozone.codes[ez]

  # Y.map.predicted.elev_p95.RF <- Y.map.predicted.elev_p95.YAI <- Y.map.predicted.percentage_first_returns_above_2m.RF <- Y.map.predicted.percentage_first_returns_above_2m.YAI <- data.table()

    # Y.val.predicted.stddev.RF <- Y.val.predicted.RF  ## initialize the same way df to store standard deviations over the RF trees (measure of uncertainty)
  
  nr.clusters <- min(detectCores(), length(params4a$mapped.years))
  cl <- makeCluster(nr.clusters)
  registerDoParallel(cl)
  melted.predictions <- foreach (yr = 1:length(params4a$mapped.years), .combine='rbind', .multicombine=T, .packages=list.of.packages) %dopar% {
  # for (yr in 1:length(params4a$mapped.years)) {
    
    ## Load dataset from CSV
    file.name <- file.path(predictors.dir, sprintf("%s_%s_%s.csv", ecozone.code, ecozone.name, params4a$mapped.years[yr]), fsep = .Platform$file.sep)
    yearly.predictors <- fread(file.name, header = T) 
    yearly.predictors <- na.roughfix(yearly.predictors) ## substitute NAs with column median
    yearly.predictors[, Ch_attr:=factor(Ch_attr, levels=params3$Ch_attr.labels)]  ## set as factor with same levels as training set
    yearly.predictors[, Elev:=as.numeric(Elev)]  
    
    ## trick to be used in the testing phase only to have smaller datasets and a faster debugging
    if (params4a$subsetting) {
      set.seed(paramsGL$global.seed)
      pix.idx <- sample(1:nrow(yearly.predictors), params4a$nr.pts.subset)
      yearly.predictors <- yearly.predictors[pix.idx, ]
    }
    
    Y.map.predicted.RF <- Y.map.predicted.YAI <- data.table(rep(params4a$mapped.years[yr], nrow(yearly.predictors)))
    
#### RANDOM FOREST ---------------------------------------------------------------------
    
    if ("RF" %in% params4a$methods) {  ## run only if RF is in params4a$methods
       
      for (targ in params4a$targ.names.temporal.lg) {
      
          RF.model.path <- file.path(models.subdir, sprintf("RF_%s.Rdata", targ), fsep=.Platform$file.sep)
          load(RF.model.path)
          
          ## prediction on map pixels
          prediction.res <- predict(rf.RF, yearly.predictors[, stats3$final.predictors, with=FALSE], type="response", predict.all=T, nodes=F)
          
          Y.map.predicted.RF <- cbind(Y.map.predicted.RF, prediction.res$aggregate)
          
          # cmd <- sprintf("Y.map.predicted.%s.RF <- prediction.res$aggregate", targ)  ## save actual predictions in the corresponding column for each resp. variable
          # eval(parse(text=cmd))
          # Y.val.predicted.stddev.RF[,targ] <- apply(prediction.res$individual, 1, sd)  ## save std dev over the RF trees (std dev of the individual prediction of each tree, a row of prediction.res$individual)
          
      }  ## end for on params4a$targ.names.temporal.lg
      
      setnames(Y.map.predicted.RF, c("Year", params4a$targ.names.temporal.lg))
      
    }  ## end if on RF
        
#### YAIMPUTE -------------------------------------------------------------------
        
    if ("YAI" %in% params4a$methods) {  ## run only if "YAI" is specified in methods
      
      ids.val.predicted <- as.data.table( newtargets(yai.rf, yearly.predictors[, stats3$final.predictors, with=FALSE], k=1)["neiIdsTrgs"] )
      ids.val.predicted[, idx:=seq(from=1, to=nrow(yearly.predictors))]
      colnames(ids.val.predicted)[1] <- "Predicted_FCID"
      Y.map.imputed <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
      setkey(Y.map.imputed, idx)
      Y.map.imputed[, c("Predicted_FCID", "idx"):=NULL]
      Y.map.predicted.YAI <- cbind(Y.map.predicted.YAI, Y.map.imputed)
      
    } ## end if on YAI
  
    setnames(Y.map.predicted.YAI, c("Year", unlist(params3$targ.names.lg)))
    
    list(Y.map.predicted.RF, Y.map.predicted.YAI)
    
  }  ## end for on params4a$mapped.years
  
  stopCluster(cl)
  
  melted.predictions.RF <- rbindlist(melted.predictions[,1])
  melted.predictions.YAI <- rbindlist(melted.predictions[,2])
  
  file.name.melted.predictions.RF <- file.path(predicted.values.dir, sprintf("persistant_treed_RF_%s.csv", params4a$mapped.ecozones[ez]), fsep = .Platform$file.sep)
  file.name.melted.predictions.YAI <- file.path(predicted.values.dir, sprintf("persistant_treed_YAI_%s.csv", params4a$mapped.ecozones[ez]), fsep = .Platform$file.sep)
  
  fwrite(melted.predictions.RF, file.name.melted.predictions.RF)
  fwrite(melted.predictions.YAI, file.name.melted.predictions.YAI)
  
}  ## end for on params4a$mapped.ecozones


#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog4a, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
