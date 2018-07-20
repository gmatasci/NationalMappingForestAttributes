## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca)
## File Name: prog6a_predict_on_random_pix.R                           
## Objective: Predict forest attributes through time for randomly sampled pixels whose predictors are provided in a CSV file.

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - set appropriate FINAL_MODEL folder

# - currently, in its parallel version, it works only with 1 method, YAI
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

param_file <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

base_wkg_dir <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/wkg"

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

param_file_prog3 = file.path(base_wkg_dir, 'AllUTMzones_params3.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog3)

source("D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/code/Functions_NatMapping_R.R")

#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params6a <- list()

params6a$experiment.name <- "NO_CHATTR"

params6a$subsetting <- T  ## to subset the dataset by a factor = params6a$data.subs.factor (for debugging in development phase or to work with smaller datasets)
# params6a$data.subs.factor <- 3
# params6a$data.subs.factor <- 1.5
params6a$data.subs.factor <- 300
params6a$nr.clusters <- 1   ## only works with 1 cluster and %do% in foreach()

# params6a$temporal.base.dir <- "E:/NTEMS/UTM_results/2017_02_27_SAMPLES_treed_changes_plots"  ## subdirectory to save files with predicted values
params6a$temporal.base.dir <- "D:/Postdoc_Giona_2016-2017/DATA_E/NTEMS/UTM_results/samples_1984_2016"

# params6a$temporal.data.types <- c("raw", "fitted")   ## type of temporal data to be analyzed: "raw" is the dataset with raw BAP proxy values, "fitted" is the dataset with fitted spectral trends
params6a$temporal.data.types <- c("fitted")

# params6a$change.types <-  c("fire", "harvesting", "always_treed")
# params6a$change.types <-  c("fire", "harvesting")
params6a$change.types <-  c("always_treed")
# params6a$change.types <-  c("harvesting")

params6a$models.subdir <- "D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/FINAL_RUN_fitted_NO_CHATTR/Models"  ## model subdirectory to load models from

params6a$mapped.years <- seq(from=1984, to=2016)
# params6a$mapped.years <- seq(from=1984, to=1986)

params6a$methods <- list("YAI")
# params6a$methods <- list("RF", "YAI")   ## has to be set as list("RF", "YAI"), that is both methods have to be run

params6a$mapped.ecozones <- list("Atlantic Maritime"=2, "Boreal Cordillera"=3, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Hudson Plains"=7, "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Cordillera"=15,
                                 "Taiga Plains"=16, "Taiga Shield East"=17, "Taiga Shield West"=18)

# params6a$harvest.ecozones <- list("Atlantic Maritime"=2, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
#                                  "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Plains"=16)

params6a$harvest.ecozones <- list("Boreal Shield East"=5, "Boreal Shield West"=6,
                                  "Montane Cordillera"=9, "Pacific Maritime"=11, "Taiga Plains"=16)


params6a$fire.ecozones <- list("Boreal Cordillera"=3, "Boreal Plains"=4, "Boreal Shield East"=5, "Boreal Shield West"=6,
                                 "Hudson Plains"=7, "Montane Cordillera"=9, "Taiga Cordillera"=15,
                                 "Taiga Plains"=16, "Taiga Shield East"=17, "Taiga Shield West"=18)

# params6a$harvest.ecozones <- list("Montane Cordillera"=9)
# 
# params6a$fire.ecozones <- list("Montane Cordillera"=9)

params6a$targ.names.temporal.lg <- c("elev_p95", "percentage_first_returns_above_2m", "gross_stem_volume", "total_biomass", "loreys_height", "basal_area")

param_file_prog6a = file.path(base_wkg_dir, 'AllUTMzones_params6a.Rdata', fsep = .Platform$file.sep)
save(params6a, file = param_file_prog6a)

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("xtable",   ## to print Latex tables
                      "sp",
                      "spdep",
                      "ggplot2",
                      "gridExtra",
                      "psych",
                      "plyr",
                      # "dplyr",   ## to be loaded before foreach to avoid "assertion failed" errors
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
if(length(new.packages)){ 
   install.packages(new.packages)
}
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}


#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time


#### LOAD MODELS & READ DATA ------------------------------------------------

pix.idx.list <- list()  ## initialize list containing subsampling indices (varying only by change.type and ecozone)
for (temporal.data.type in params6a$temporal.data.types) {
  
  ## Load model files either in the "raw" or "fitted" folder
  load(file.path(params6a$models.subdir, sprintf('%s_%s', temporal.data.type, params6a$experiment.name), "YAI.Rdata", fsep=.Platform$file.sep))
  load(file.path(params6a$models.subdir, sprintf('%s_%s', temporal.data.type, params6a$experiment.name), "Ytrn_longID.Rdata", fsep=.Platform$file.sep))
  
  stats_file_prog3 <- sprintf("D:/Postdoc_Giona_2016-2017/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR/FINAL_RUN_fitted_NO_CHATTR/results_%s_%s/stats3.Rdata", temporal.data.type, params6a$experiment.name)
  load(stats_file_prog3)
  
  print(sprintf("Temporal data: %s", temporal.data.type))

  for (change.type in params6a$change.types) {
    
    print(sprintf("Change type: %s", change.type))
    
    samples.dir <- file.path(params6a$temporal.base.dir, change.type, 'samples_u', fsep = .Platform$file.sep) 
    
    if (temporal.data.type == "fitted") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics_fitted', fsep = .Platform$file.sep)
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, sprintf('predicted_values_fitted_%s', params6a$experiment.name), fsep = .Platform$file.sep)
    } else if (temporal.data.type == "raw") {
      predictors.dir <- file.path(params6a$temporal.base.dir, change.type, 'metrics', fsep = .Platform$file.sep) 
      predicted.values.dir <- file.path(params6a$temporal.base.dir, change.type, sprintf('predicted_values_%s', params6a$experiment.name), fsep = .Platform$file.sep) 
    }
    
    if (! file.exists(predicted.values.dir)){dir.create(predicted.values.dir, showWarnings = F, recursive = T)}
    
    ## Select the appropriate list of ecozones wrt ch_attr
    if (change.type == "always_treed") {
      ecozones <- params6a$mapped.ecozones
    } else if (change.type == "fire") {
      ecozones <- params6a$fire.ecozones
    } else if (change.type == "harvesting") {
      ecozones <- params6a$harvest.ecozones
    }
    
    nr.zero.samples.removed <- data.frame(matrix(nrow=length(ecozones), ncol=length(params6a$mapped.years)))
    rownames(nr.zero.samples.removed) <- names(ecozones)
    colnames(nr.zero.samples.removed) <- paste0('yr', params6a$mapped.years)
    pct.aberr.elev <- nr.zero.samples.removed
    for (ez in 1:length(ecozones)) {
      
      ecozone.name <- names(ecozones[ez])
      ecozone.code <- ecozones[[ez]]
      
      print(ecozone.name)
      
      file.name <- file.path(samples.dir, sprintf("%s_%s.csv", ecozone.code, ecozone.name), fsep = .Platform$file.sep)
      if (! file.exists(file.name)) {
        next
      }
      pix.IDs.ecoz <- fread(file.name, header=F, verbose=F)
      
      #### TO USE WITH FOR ON params6a$mapped.years ------------
      predictions.list <- lapply(1:length(params6a$methods), function(x) NULL)
      names(predictions.list) <- params6a$methods
      #### TO USE WITH FOR ON params6a$mapped.years ------------
      
      #### FOREACH ON params6a$mapped.years: not working because heavy model object has to be passed to each core, thus saturating RAM 
      # nr.clusters <- round(length(params6a$mapped.years)/8)+1
      # cl <- makeCluster(nr.clusters)
      # registerDoParallel(cl)
      # melted.preds.dt <- foreach(yr = 1:length(params6a$mapped.years), .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=T for more info when debugging
      #### FOREACH ON params6a$mapped.years
      for (yr in 1:length(params6a$mapped.years)) {

        if (yr %% 5 == 0) {
          print(params6a$mapped.years[yr])
        }
        
        ## Load dataset from CSV
        if (temporal.data.type == "fitted") {
          file.name <- file.path(predictors.dir, sprintf("%s_%s_%s_fitted.csv", ecozone.code, ecozone.name, params6a$mapped.years[yr]), fsep = .Platform$file.sep)
        } else {
          file.name <- file.path(predictors.dir, sprintf("%s_%s_%s.csv", ecozone.code, ecozone.name, params6a$mapped.years[yr]), fsep = .Platform$file.sep)
        }

        yearly.predictors <- fread(file.name, header=T, verbose=F)
        
        ## Subsample the dataset by a factor params6a$data.subs.factor
        if (params6a$subsetting) {   ## if subsampling is specified...
          str.chType <- sprintf('%s_%s', change.type, ecozone.code)  ## element name to store/retrieve subsampling pixel indices: they only vary wrt change type and ecozone
          if (yr == 1 & temporal.data.type == params6a$temporal.data.types[1]) {  ## generate random indices only for 1st year and for raw data, as all years and raw/fitted will refer to the same pixels
            set.seed(paramsGL$global.seed)
            pix.idx <- sample(1:nrow(yearly.predictors), round(nrow(yearly.predictors)/params6a$data.subs.factor))
            pix.idx.list[[str.chType]] <- pix.idx   ## save list of inidices in a dynamic list
          } else {   ## for the subsequent years or other data types
            pix.idx <- pix.idx.list[[str.chType]]   ## retrieve indices from dynamic list
          }
        } else {  ## ...if not keep all the rows
          pix.idx <- 1:nrow(yearly.predictors)
        }
        yearly.predictors <- yearly.predictors[pix.idx, ]
        pix.IDs <- pix.IDs.ecoz[pix.idx, ]
        
        ## Change negative YrsSince_GrCh values back to positive
        yearly.predictors$YrsSince_GrCh <- abs(yearly.predictors$YrsSince_GrCh)
        
        pct.aberr.elev[ez, yr] <- sum(yearly.predictors$Elev > 5000) / nrow(yearly.predictors)
     
        ## Set unique pixel ID, with UTM zone info and index value
        yearly.predictors[, pixID:=sprintf('%s_%d', pix.IDs[[1]], pix.IDs[[2]])]
        yearly.predictors[, index:=NULL]
        
        ## Preprocess columns
        yearly.predictors[, seq(1,ncol(yearly.predictors)-1):= na.roughfix(yearly.predictors[, !colnames(yearly.predictors) == 'pixID', with=FALSE])] ## substitute NAs with column median
        yearly.predictors[, Ch_attr:=factor(Ch_attr, levels=params3$Ch_attr.labels)]  ## set as factor with same levels as training set
        yearly.predictors[, Elev:=as.numeric(Elev)]

        samples.before <- nrow(yearly.predictors)
        yearly.predictors <- yearly.predictors[b1!=0 & b2!=0 & b3!=0 & b4!=0 & b5!=0 & b6!=0]
        nr.zero.samples.removed[ez, yr] <- samples.before - nrow(yearly.predictors)

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
          
          #### FOREACH ON 1 SPLIT OF yearly.predictors: temporary version working only with %do% and nr.clusters = 1
          nr.rows <- nrow(yearly.predictors)
          Y.map.predicted.YAI <- data.table(rep(params6a$mapped.years[yr], nr.rows))
          nr.clusters <- params6a$nr.clusters
          split.vect <- rep(1, nr.rows)   ## split dataframe into one chunk only
          yearly.predictors.splits <- split(yearly.predictors, split.vect)
          cl <- makeCluster(nr.clusters)
          registerDoParallel(cl)
          ## CHANGE TO %dopar% FOR PARALLEL (NOT WORKING)
          Y.map.imputed <- foreach(split = 1:nr.clusters, .combine='rbind', .packages=list.of.packages) %do% {   #add .verbose=T for more info when debugging
            ids.val.predicted <- as.data.table( newtargets(yai.rf, yearly.predictors.splits[[split]][, stats3$final.predictors, with=FALSE], k=1)["neiIdsTrgs"] )
            ids.val.predicted[, order:=seq(from=1, to=nrow(yearly.predictors.splits[[split]][, stats3$final.predictors, with=FALSE]))]   ## add index telling the order of the samples before the merge
            colnames(ids.val.predicted)[1] <- "Predicted_FCID"
            Y.map.imputed.split <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
            setkey(Y.map.imputed.split, order)   ## after the merge resort with respect to that index...
          }
          stopCluster(cl)
          
          Y.map.imputed[, c("Predicted_FCID", "order"):=NULL]   ## ...and then remove it

          Y.map.predicted.YAI <- cbind(Y.map.predicted.YAI, Y.map.imputed)

          setnames(Y.map.predicted.YAI, c("Year", unlist(params3$targ.names.lg)))
          
          #### TO USE WITH FOREACH ON params6a$mapped.years ------------
          # cbind(pixID=yearly.predictors[,pixID], Y.map.predicted.YAI)
          #### TO USE WITH FOREACH ON params6a$mapped.years ------------
          
          #### TO USE WITH FOR ON params6a$mapped.years ------------
          predictions.list$YAI <- rbind(predictions.list$YAI, list(Y.map.predicted.YAI))
          #### TO USE WITH FOR ON params6a$mapped.years ------------
          
        } ## end if on YAI
      
      }  ## end for on params6a$mapped.years
      
      #### TO USE WITH FOREACH ON params6a$mapped.years ------------
      # stopCluster(cl)
      # melted.predictions <- list()
      # melted.predictions$YAI <- melted.preds.dt
      #### TO USE WITH FOREACH ON params6a$mapped.years ------------
      
      #### TO USE WITH FOR ON params6a$mapped.years ------------
      ## Stitch together the list elements by column of lists
      melted.predictions <- lapply(predictions.list, function(x) rbindlist(x))
      
      rm(predictions.list)
      gc()

      ## Add pixel ID replicated by year to later reshape the data with one row per pixel
      melted.predictions <- lapply(melted.predictions, cbind, pixID=rep(yearly.predictors[,pixID], length(params6a$mapped.years)))
      #### TO USE WITH FOR ON params6a$mapped.years ------------
      
      ## Apply reshape_save() function on every column of melted.predictions to obtain a dt with pixels in rows and years in columns (each row is a pixel attribute time-series)
      lapply(seq_along(melted.predictions), function(idx) reshape_save(melted.predictions[[idx]], predicted.values.dir, change.type, names(melted.predictions)[idx], params6a$targ.names.temporal.lg, params6a$mapped.years, ecozone.code, ecozone.name))
      
      rm(melted.predictions)
      gc()
      
    }  ## end for on ecozones
    
    ## ---------- UNCOMMENT ------------
    zeroes.file <- file.path(params6a$temporal.base.dir, sprintf('ZeroSamplesRemovedByEcoByYear_%s.Rdata', temporal.data.type), fsep = .Platform$file.sep) 
    save(nr.zero.samples.removed, file = zeroes.file)
    ## ---------- UNCOMMENT ------------
    
    aberr.elev.file <- file.path(params6a$temporal.base.dir, sprintf('AberrElevByEcoByYear_%s.csv', temporal.data.type), fsep = .Platform$file.sep) 
    write.csv(pct.aberr.elev, aberr.elev.file, row.names=T)
  
  }  ## end for on params6a$change.types
  
}  ## end for on params6a$temporal.data.types


#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog6a, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)


