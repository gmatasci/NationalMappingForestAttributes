#### CODE INFOS -------------------------------------------------------------------

## Project Name: NationalImputationForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name:                            
## Objective: Descriptive statistics and plots, variable selection and Random Forest prediction/imputation

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - check final UTM zones to sample (now 17 with the removal of 11S)
# - params3a$parallel.RF <- T
# - params3a$subsetting <- F

# - run with final data: comparison between predicting response variables one at a time and not as a single multivariate Y: 
#   this means comparing yaimpute rf, rf, etc.
# - ############# = TODEL
# - to check to see sampled values based on "FCID" or "POLY250ID" from X or Y dataframes loaded here linked to field "POLY250ID" of <UTMzone>_pt_centerpt_training_validation.shp (saved in /BAP_Imputation_working/wkg/<UTMzone>/)
#   then related to <UTMzone>_plot_inventory_attributes2.csv by "unique_id"

## SOLVED:
# -V write script between prog2 and prog3 to group data by ecozone rather than by UTM zone -- no, we run the analysis at the national level with one single model (long, lat & other trends should account)
# -V change loop over ecozones instead of over UTM zones -- no need to do it as we use one single model valid for all Can territory
# -V multi.hist() per var by UTM zone -- no, only works for dataframes with same nr of columns (UTM zones have different nr of samples)
# -V ANOVA style qplot for categorical vars -- boxplot of key targets by category for CAN and barplot by UTM zone
# -V ggpair each type vs the targets: fix loop for dynamic plots (consider copying the same 4 times) -- done with pairs.panels()
# -V check plots section with new structure
# -V test density=T for pair.panels -- no, this just adds a density curve to the histogram in the diagonal
# -V global correlation matrix without subsampling to see redundant variables -- added with filters for high and low correlations
# -V subsample only for pair.plot and not for hist or rest -- done!
# -V consider setting mtry=sqrt(nr var) -- set as parameter to be used
# -V 10 run and report mean+-std?  -- not done by harold and since it is not a model comparison paper we can drop it (minimal variation around the mean) and just use same seed each time
# -V implement model selection based on Impute_Vs_predict.R -- done
# -V build model with proprer data -- final statistical dataset built on 19/05/2016
# -V make sure Ch_attr is interpreted as a categorical/factor variable by RF and that all the factor values are found in TRN (no new values in VAL), otherwise errors (http://stats.stackexchange.com/questions/29446/random-forest-and-new-factor-levels-in-test-set)
# -V yaImpute uses a RF with a nrTrees shared across all Ys, so actual nrTrees is nrTrees/nrYs -- solved by setting ntree=params3a$ntree*length(params3a$targ.yaImp.names.lg) in YaImp parameters
# -V remove ECOZONE NA filtering once prog0 to 2 are re-run
# -V uncomment lines related to Ch_attr
# -V check descriptive statistics part with full data


#### READS/WRITES ------------------------------------------------------------

## READS:
# - "lidar_metrics_mean_training_validation.csv": (from prog1) csv table with all observed LiDAR and forest attributes (Y) for selected TRN and VAL samples (3x3 polygons with average plot values or just single plot value) for all UTM zones
# - "poly_training_validation_exvars_extract.csv":  (from prog2) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

## WRITES:
# - "yai_rf_gnn_sk.RData": R data file with yaimpute models and training set to use later for mapping (prediction/imputation) on the grid

#### INIT --------------------------------------------------------------------

print('Prog3a: descriptive stats/plots, model selection and model assessment') 

rm(list=ls())

param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

source("D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/Functions_NatImp.R")

#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params3a <- list()

## final target variables to be predicted by both RF and YaImp
params3a$targ.names.lg <- list("elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean", "loreys_height", "basal_area", "gross_stem_volume", "total_biomass")
params3a$targ.names.sh <- list("el_m","el_std","el_p95","el_cv","pct_1r_ab2","pct_1r_abm", "loreys_h", "basal_a", "stem_v", "tot_biom")

## target variables given as input to YaImp to find best neighbor (other 4 are attached once FCIDs are known) 
params3a$targ.yaImp.names.lg <- list("elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean")
params3a$targ.yaImp.names.sh <- list("el_m","el_std","el_p95","el_cv","pct_1r_ab2","pct_1r_abm")

## two key target variables used in descriptive plots to show X to Y relationship
params3a$targ.key.names.lg <- list("percentage_first_returns_above_2m", "total_biomass")
params3a$targ.key.names.sh <- list("pct_1r_ab2", "tot_biom")

params3a$run.checks <- F   ## where to run checks for duplicate FCIDs, run only once, as of 01/06/2016 all is OK
params3a$descriptive.stats <- F  ## wheter to run or not descriptive stats block

## retained 5 classes: 0="No change", 1="Fire", 2="Harvesting", 3="Non-stand replacing", 4="Infrastructure"
params3a$Ch_attr.classes <- c("No change", "Fire", "Harvesting", "Non-stand replacing", "Infrastructure")
params3a$Ch_attr.labels <- c(0, 1, 2, 3, 4)

params3a$groups.to.plot <- list('Bands', 'TCcomps_VI', 'ChangeAttribution', 'Years_Topo_Coords')  ## to plot relation with key response vars separately (and in a visually pleasing way)

params3a$nr.pts.plot <- 2000 ## to make pairs plots less dense in points and to subset the dataset for tests in the developement phase

params3a$min.year.of <- 1960  ## used to replace the 0 in year.of (no change areas are set as last changed in this year) and produce meaningful scatterplots (the lower the year the more developed the canopy)

params3a$corr.high <- 0.9  ## corr threshold to highlight high correlation among predictors
params3a$corr.low <- 1  ## corr threshold to highlight low correlation between predictors and response vars

params3a$nr.of.bins <- 12  ## nr of bins for UTM zone histograms
params3a$hist.centile <- 0.999  ## percentile threshold to subset data for histograms (avoid extreme values setting a too high xlim)

## TO CHANGE IF EXPERT-BASED APPROACH IS USED
# params3a$redund.pred.names.sh <- list("b4", "b5", "TCD", "TCA", "EVI")  ## variables selected for removal after inspection of the Corr matrix bc redundant
## TO CHANGE IF EXPERT-BASED APPROACH IS USED

params3a$predictor.groups <- list("all", "SKstudy", "Rthresh_0p8", "Rthresh_0p95", "Rthresh_0p9_no_coords", "Rthresh_0p9_no_change", "Rthresh_0p9")
params3a$best.predictor.group <- "Rthresh_0p9"

params3a$all.ecozones <- list("Arctic Cordillera", "Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", 
                          "Hudson Plains", "Mixedwood Plains", "Montane Cordillera", "Northern Arctic", "Pacific Maritime", "Semiarid Prairies", 
                          "Southern Arctic", "Subhumid Prairies", "Taiga Cordillera", "Taiga Plains", "Taiga Shield East", "Taiga Shield West")

params3a$sampled.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", "Hudson Plains", 
                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")  ## Semiarid Prairies and Subhumid Prairies (ECOZONE_ID = 10 in shp) result in 0 and 3 plots, respectively, so we remove them

## Actual parameters to be used
params3a$run.MS <- T
params3a$run.MA <- T
params3a$parallel.RF <- T
params3a$subsetting <- F
params3a$ntree.MS <- 32   ## ideally 500 but then YaImp (not parallel) takes ages
params3a$ntree.MA <- 500  

## Parameters for testing
# params3a$run.MS <- T
# params3a$run.MA <- T
# params3a$parallel.RF <- F
# params3a$subsetting <- T
# params3a$ntree.MS <- 20
# params3a$ntree.MA <- 50

params3a$mtry <- 'sqrt_nr_var'
params3a$nodesize <- 5

params3a$metrics.MS <- "rsq"
params3a$metrics.MA <- c("ymean", "yrange", "ystdev", "xmean", "xrange", "xstdev", "rsq", "rmsd", "nrmsd", "ac_uns", "ac_sys", "bias")
params3a$metrics.MA.colnames <- c("o_mean", "o_rangemin", "o_rangemax", "o_stdev", "p_mean", "p_rangemin", "p_rangemax", "p_stdev", "rsq", "rmsd", "nrmsd", "ac_uns", "ac_sys", "bias")


param_file_prog3a = file.path(base_wkg_dir, 'AllUTMzones_params3a.Rdata', fsep = .Platform$file.sep) 
save(params3a, file = param_file_prog3a)

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("rgdal",
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
                      "data.table"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time

#### READ DATA --------------------------------------------------------------

## load plot training and validation data for lidar metrics and derived structural attributes (old "str", for structural)
Y.trn.val.raw <- read.csv(file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", fsep = .Platform$file.sep))

## change units for "summable" variables: kg per plot (625m^2) to kg per hectare (10000m^2) --> multiply by 16
Y.trn.val.raw[, c("basal_area", "gross_stem_volume", "total_biomass")] <- Y.trn.val.raw[, c("basal_area", "gross_stem_volume", "total_biomass")]*16

## load explanatory variables extracted for training plots (old "env", for environment)
X.trn.val.raw <- read.csv(file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", fsep = .Platform$file.sep))

X.trn.val.raw$Ch_attr <- as.factor(X.trn.val.raw$Ch_attr)   ## to let the RF interpret Ch_attr as a categorical variable

if (params3a$run.checks) {
  
  ## check for order of sample points in the two datasets X and Y
  if ( !identical(X.trn.val.raw$FCID, Y.trn.val.raw$FCID) ) {
    stop("Prog3a: order of FCIDs is different in lidar_metrics_mean_training_validation.csv wrt poly_training_validation_exvars_extract.csv")
  } 
  
  ## check for duplicate FCIDs in the dataset
  if ( length(unique(X.trn.val.raw$FCID)) != length(X.trn.val.raw$FCID) ) {
    stop("Prog3a: FCIDs are not unique in lidar_metrics_mean_training_validation.csv and poly_training_validation_exvars_extract.csv")
  }
  
  ## check for duplicate FCIDs (POLY250ID) in the various <UTMzone>_poly_training_validation.shp shapefiles
  uniq.df <- data.frame(matrix(ncol = 2, nrow = 0))
  FCIDs.df <- vector(length=0)
  for (z in 1:length(paramsGL$zones)) {
    zone <- paramsGL$zones[z]
    zone.nr <- substr(zone, 4, nchar(zone))
    wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)
    poly.training.validation <-  readOGR(dsn = wkg_dir, layer = paste(zone,"_poly_training_validation",sep = ''))
    uniq.df <- rbind( uniq.df, c( length(poly.training.validation@data$POLY250ID), length(unique(poly.training.validation@data$POLY250ID)) ) )
    FCIDs.df <- c(FCIDs.df, as.vector(poly.training.validation@data$POLY250ID) )
  }
  if ( colSums(uniq.df)[1] != colSums(uniq.df)[2] ) {
    stop("Prog3a: There are duplicate FCIDs (POLY250ID) within the various <UTMzone>_poly_training_validation.shp shapefiles")
  }
  if ( length(unique(FCIDs.df)) != length(FCIDs.df) ) {
    stop("Prog3a: There are duplicate FCIDs (POLY250ID) between the various <UTMzone>_poly_training_validation.shp shapefiles")
  }

}

## filter both dataframes to remove plots with 0 values in Landsat (SRef_<UTMzone>_2010_proxy_v2.dat) or Topo (<UTMzone>_DEM.dat) data due to small overlap mismatch between sampling polygons shp and raster layers
idx.Landsat.zero <- X.trn.val.raw[, "b1"]==0 & X.trn.val.raw[, "b2"]==0 & X.trn.val.raw[, "b3"]==0 & X.trn.val.raw[, "b4"]==0 & X.trn.val.raw[, "b5"]==0 & X.trn.val.raw[, "b7"]==0
idx.Elev.zero <- X.trn.val.raw[, "Elev"]==0
idx.to.remove <- idx.Landsat.zero | idx.Elev.zero  ## set as TRUE if either of the 2 conditions is TRUE
Y.trn.val <- Y.trn.val.raw[!idx.to.remove, ]  ## keep rows that are not (!) to remove
X.trn.val <- X.trn.val.raw[!idx.to.remove, ]
rm(Y.trn.val.raw, X.trn.val.raw)  ## clear space after filtering

## delete samples from Subhumid Prairies (just 3 plots after filtering)
idx.Subhumid.Prairies <- Y.trn.val$ECOZONE == "Subhumid Prairies"
Y.trn.val <- Y.trn.val[!idx.Subhumid.Prairies, ]
X.trn.val <- X.trn.val[!idx.Subhumid.Prairies, ]
Y.trn.val$ECOZONE <- factor(Y.trn.val$ECOZONE)  ## to drop the level we just removed

## save stats about removed samples and nr of sample from each Ecozone, UTM zone and Ch_attr
stats3a <- list()
stats3a$NA.samples.removed <- sum(idx.to.remove==T) + sum(idx.Subhumid.Prairies==T)
TRN.samples.per.Ecozone <- table(Y.trn.val[Y.trn.val$TV=="TRAINING",]$ECOZONE)  ## to be added later on to the final stats dataframes
VAL.samples.per.Ecozone <- table(Y.trn.val[Y.trn.val$TV=="VALIDATION",]$ECOZONE)
TRN.samples.per.UTMzone <- table(X.trn.val[X.trn.val$TV=="TRAINING",]$UTMzone)
VAL.samples.per.UTMzone <- table(X.trn.val[X.trn.val$TV=="VALIDATION",]$UTMzone)
TRN.samples.per.Ch_attr <- table(X.trn.val[X.trn.val$TV=="TRAINING",]$Ch_attr)
VAL.samples.per.Ch_attr <- table(X.trn.val[X.trn.val$TV=="VALIDATION",]$Ch_attr)

#### CORRELATION MATRIX ----------------------------------------------

## compute the following outside of if block bc we want to know redundant predictors in all cases
cont.idx <- !unlist(params2$pred.names.sh) %in% c(params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr) ## subset to continuous predictors
X.trn.val.corr <- X.trn.val[, unlist(params2$pred.names.sh)[cont.idx]]
X.corr.matrix <- cor(X.trn.val.corr)  ## compute corr only among predictors (top-left sub-block)

## find subset of redundant variables automatically: in a pair having r>cutoff remove the variable with highest average absolute correlation with rest of available variables (exact=T)
stats3a$redund.pred.names.p80 <- findCorrelation(X.corr.matrix[seq(1,nrow(X.corr.matrix)-2), seq(1,nrow(X.corr.matrix)-2)], cutoff=0.80, exact=T, verbose=F, names=T)  ## do not consider 2 last variables to avoid excluding Lat/Long
stats3a$redund.pred.names.p90 <- findCorrelation(X.corr.matrix[seq(1,nrow(X.corr.matrix)-2), seq(1,nrow(X.corr.matrix)-2)], cutoff=0.90, exact=T, verbose=F, names=T)
stats3a$redund.pred.names.p95 <- findCorrelation(X.corr.matrix[seq(1,nrow(X.corr.matrix)-2), seq(1,nrow(X.corr.matrix)-2)], cutoff=0.95, exact=T, verbose=F, names=T)

if (params3a$descriptive.stats) {
  
  Y.trn.val.corr <- Y.trn.val[, unlist(params3a$targ.names.lg)]
  colnames(Y.trn.val.corr) <- params3a$targ.names.sh
  
  X.corr.matrix.high <- X.corr.matrix
  X.corr.matrix.high[abs(X.corr.matrix) < params3a$corr.high] <- NA  ## keep only high correlations among predictors
  
  data.matrix <- cbind(X.trn.val.corr, Y.trn.val.corr)
  corr.matrix <- cor(data.matrix)
  XtoY.corr.matrix <- corr.matrix[1:ncol(X.trn.val.corr), seq(ncol(X.trn.val.corr)+1,ncol(corr.matrix),1)]  ## subset to keep only top-right corr matrix sub-block  
  XtoY.corr.matrix.low <- XtoY.corr.matrix
  XtoY.corr.matrix.low[abs(XtoY.corr.matrix) > params3a$corr.low] <- NA   ## keep only low predictor-response correlations
  
  XY.corr.matrix <- round(cbind(X.corr.matrix, XtoY.corr.matrix), digits=3)  ## put together the two matrix sub-blocks
  XY.corr.matrix.NA <- round(cbind(X.corr.matrix.high, XtoY.corr.matrix.low), digits=3)
  
  write.csv(XY.corr.matrix, file = file.path(base_results_dir, "XY_corr_matrix.csv", sep = ''))
  str <- sprintf( "XY_corr_matrix_filtr_High%s_Low%s.csv", gsub(".", "p", as.character(params3a$corr.high), fixed=T), gsub(".", "p", as.character(params3a$corr.low), fixed=T)) 
  write.csv(XY.corr.matrix.NA, file = file.path(base_results_dir, str, sep = ''))

#### CANADIAN LEVEL X-Y RELATION PLOTS ------------------------------------------
  
  ## create figure subdirectories
  CAN.subdir <- file.path(base_figures_dir, "Relations_CAN_level", fsep = .Platform$file.sep)
  UTMzone.subdir <- file.path(base_figures_dir, "Histograms_UTMzone_level", fsep = .Platform$file.sep)
  if (! file.exists(CAN.subdir)){dir.create(CAN.subdir, showWarnings = F, recursive = T)}
  if (! file.exists(UTMzone.subdir)){dir.create(UTMzone.subdir, showWarnings = F, recursive = T)}
  
  ## create indices for subsampling 
  set.seed(paramsGL$global.seed)
  
  plot.idx.init <- sample(1:nrow(X.trn.val), params3a$nr.pts.plot)  # take a subset of # params3a$nr.pts.plot of points to plot
  UTM.idx <- which(paste('UTM', X.trn.val$UTMzone, sep='') %in% paramsGL$zones)  ## to only select samples of the UTM zones specified in paramsGL$zones
  plot.idx <- plot.idx.init[plot.idx.init %in% UTM.idx]
  
  Y.data.to.plot <- Y.trn.val[plot.idx, unlist(params3a$targ.key.names.lg)] ## select Y columns to plot
  colnames(Y.data.to.plot) <- params3a$targ.key.names.sh
  
  gr.idx <- 1
  for (pred.gr.names in params3a$groups.to.plot) {
    gr.name <- unlist(params3a$groups.to.plot[gr.idx])
    gr.idx <- gr.idx+1
    
    if (gr.name == 'ChangeAttribution') {
      pred.names <- "Ch_attr"
      X.data.to.plot <- as.data.frame(X.trn.val[, pred.names])
      colnames(X.data.to.plot) <- pred.names
      
      Y.data.to.plot.Ch_attr <- Y.trn.val[, unlist(params3a$targ.key.names.lg)]
      colnames(Y.data.to.plot.Ch_attr) <- params3a$targ.key.names.sh
      
      df.to.plot <- cbind(X.data.to.plot, Y.data.to.plot.Ch_attr)

      ## Canadian level (merged UTM zones) boxplots by Change_attr (categorical var)
      plot1 <- ggplot(df.to.plot, aes(x=Ch_attr, y=pct_1r_ab2, fill=Ch_attr)) + 
        geom_boxplot(notch=F, outlier.shape = NA) + 
        scale_x_discrete(labels=params3a$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      plot2 <- ggplot(df.to.plot, aes(x=Ch_attr, y=tot_biom, fill=Ch_attr)) + 
        geom_boxplot(notch=F, outlier.shape = NA) + 
        scale_x_discrete(labels=params3a$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylim(0, 20000)
      str <- file.path(CAN.subdir, sprintf("CAN_boxplots_keyTarg_VS_%s.pdf", gr.name), sep='')
      pdf(str, width=10, height=5)
        grid.arrange(plot1, plot2, nrow=1, ncol=2, top=sprintf("%s VS pct_1ret_ab_2m and tot_biom", gr.name))
      dev.off()
      
      next
      
    } else if (gr.name == 'Bands') {
      pred.names <- unlist(params2$pred.names.sh$bands)
      cex.val <- 1
      cex.labels.val <- 1
      cex.cor.val <- 0.8
    } else if (gr.name == 'TCcomps_VI') {
      pred.names <- c(unlist(params2$pred.names.sh$TC), unlist(params2$pred.names.sh$VI))
      cex.val <- 1
      cex.labels.val <- 1
      cex.cor.val <- 0.8
    } else if (gr.name == 'Years_Topo_Coords') {
      pred.names <- c(unlist(params2$pred.names.sh$yrs), unlist(params2$pred.names.sh$topo), unlist(params2$pred.names.sh$coords))
      cex.val <- 1
      cex.labels.val <- 1
      cex.cor.val <- 0.8
    }
    
    ## Canadian level (merged UTM zones) bivariate scatterplots with linear model R on it
    X.data.to.plot <- X.trn.val[plot.idx, pred.names]
    data.to.plot <- cbind(X.data.to.plot, Y.data.to.plot)
    if (gr.name == 'Years_Topo_Coords') {
      data.to.plot$GreatCh_yr[data.to.plot$GreatCh_yr == 0] <- params3a$min.year.of
    }
    str <- file.path(CAN.subdir, sprintf("CAN_pairsPanel_keyTarg_VS_%s.pdf", gr.name), sep='')
    pdf(str)
      pairs.panels(data.to.plot, pch=1, cex=cex.val, 
                   scale=F, density=F, ellipses=F, lm=T, jiggle=F, rug=F,
                   breaks=15, cex.labels=cex.labels.val, cex.cor=cex.cor.val,
                   main=sprintf("%s VS pct_1ret_ab_2m and tot_biom", gr.name))
    dev.off()
  
  }
  

#### HISTOGRAMS BY UTM ZONE ------------------------------------------
  
  ## continuous var
  cat.vars <- c("GreatCh_yr", "Ch_attr")  ## categorical vars to remove
  cont.var.idx <- !(colnames(X.trn.val) %in% cat.vars)  ## get indixes of all predictors except "Ch_attr" and "GreatCh_yr" (dealt separately as a barplot)
  data.to.plot <- cbind(X.trn.val[, cont.var.idx], Y.trn.val[, unlist(params3a$targ.names.lg)])  ## not subset bc histograms can handle many points 
  for (v in 6:ncol(data.to.plot)) {   ## from column 6 onwards (5 is UTMzone and is the last one before the predictors)
    var <- colnames(data.to.plot)[v]
    var.to.plot.compl <- data.to.plot[, v]
    quant.idx <- var.to.plot.compl < quantile(var.to.plot.compl, params3a$hist.centile, names=F) & var.to.plot.compl > quantile(var.to.plot.compl, 1-params3a$hist.centile, names=F) ## index to select data within the 1-99th centile range to avoid extreme values skewing the histogram
    var.to.plot <- var.to.plot.compl[quant.idx] 
    UTMzone.to.plot <- data.to.plot[quant.idx, 5]
    my.hist.lims <- c(min(var.to.plot), max(var.to.plot))  
    my.breaks <- c( seq(min(var.to.plot), max(var.to.plot), l=params3a$nr.of.bins+1) )
    str <- file.path(UTMzone.subdir, sprintf("UTMzonesHist_%s.pdf", var), sep='')
    pdf(str)
    par(mfrow=c(4,5), oma = c(5,4,2,2) + 0.1, mar = c(2,2,2,2) + 0.1)  ## parameters to arrange nicely the panel of plots
    for (z in 1:length(paramsGL$zones)) {
      zone <- paramsGL$zones[z]
      zone.nr <- substr(zone, 4, nchar(zone))
      zone.idx <- UTMzone.to.plot == zone.nr  ## get zone indexes to plot only data of that specific zone in this round of the loop
      hist(var.to.plot[zone.idx], main=sprintf("%s", zone), freq=T, breaks=my.breaks, xlim=my.hist.lims, xlab=NULL, ylab=NULL, col='lightgreen')  # , breaks=20
    }
    title(var, outer=T)
    dev.off() 
  }
      
  ## categorical vars
  for (var in cat.vars) {   ## from column 6 onwards (5 is UTMzone and is the last one before the predictors)
    data.to.plot <- cbind(X.trn.val[, 1:5], X.trn.val[, var])   ## not subset bc barplots can handle many points 
    str <- file.path(UTMzone.subdir, sprintf("UTMzonesHist_%s.pdf", var), sep='')
    pdf(str)
    par(mfrow=c(4,5), oma = c(5,4,2,2) + 0.1, mar = c(1,0,4.5,2) + 0.1)
    for (z in 1:length(paramsGL$zones)) {
      zone <- paramsGL$zones[z]
      zone.nr <- substr(zone, 4, nchar(zone))
      zone.idx <- which(data.to.plot$UTMzone == zone.nr)
      if (var == "GreatCh_yr") {
        yr.range <- 1985:2012   ## full range to be barplotted
        yr.padded <- c( data.to.plot[zone.idx, 6], yr.range[!yr.range %in% unique(data.to.plot[zone.idx, 6])] )  ## add to the year data the years that are missing to have a complete set 
        barplot(table(yr.padded), main=sprintf("%s", zone), ylim=c(0, 450))  ## years never appearing will have a value of 1 instead of 0 in the bar plot 
      } else {
        barplot(table(data.to.plot[zone.idx, 6]), main=sprintf("%s", zone), ylim=c(0, 1200))
      }
    }
    title(var, outer=T)
    dev.off()
  }

}

#### DATA PREPARATION -------------------------------------------------

## remove unneeded vars in X.trn.val data
X.trn.val.data <- X.trn.val[, c('FCID', unlist(params2$pred.names.sh), 'TV')]

## extract core str response variables using in random forest analysis (leave out ECOZONE, to be retrieved later on through the FCID)
## imputation map only contains single band of plot FCIDs bc extracted key variables that get attached to map after imputation
Y.trn.val.data <- Y.trn.val[, c("FCID", unlist(params3a$targ.names.lg), "TV")]

X.trn.val.data <- FCID_2_rownames(X.trn.val.data)
Y.trn.val.data <- FCID_2_rownames(Y.trn.val.data)

if (params3a$subsetting) {
  plot.idx <- sample(1:nrow(X.trn.val.data), params3a$nr.pts.plot)
  X.trn.val.data <- X.trn.val.data[plot.idx, ]
  Y.trn.val.data <- Y.trn.val.data[plot.idx, ]
}

X.trn <- subset(X.trn.val.data, TV == "TRAINING")
Y.trn <- subset(Y.trn.val.data, TV == "TRAINING")
X.val <- subset(X.trn.val.data, TV == "VALIDATION")
Y.val <- subset(Y.trn.val.data, TV == "VALIDATION")

## remove column with TV label
X.trn$TV <- NULL
Y.trn$TV <- NULL
X.val$TV <- NULL
Y.val$TV <- NULL

#### MODEL SELECTION -------------------------------------------------

if (params3a$run.MS) {
  
  print('Model selection')
  
  ## loop to choose best subset of predictors
  RES_MS.df <- data.frame(matrix(ncol=length(params3a$targ.names.lg)+1, nrow=0))
 
  case.text <- sprintf("Method (ntree=%s, mtry=%s, nodesize=%s)", params3a$ntree.MS, params3a$mtry, params3a$nodesize)
  colnames(RES_MS.df) <- c(case.text, unlist(params3a$targ.names.lg))

  for (predictor.gr in params3a$predictor.groups) {
  
    print(sprintf('Predictors group: %s', predictor.gr))
    temp.tic <- proc.time() # start clocking time for each subset

    if (predictor.gr == "all") {
      predictors <- unlist(params2$pred.names.sh)
    } else if (predictor.gr == "SKstudy") {
      predictors <- c("TCB", "TCG", "TCW", "TCA", "TCD", "GreatCh_yr", "Ch_attr", "Elev", "Slope", "TWI", "TSRI", "Long", "Lat")  ## 13 predictors
    } else if (predictor.gr == "Rthresh_0p8") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% stats3a$redund.pred.names.p80]
    } else if (predictor.gr == "Rthresh_0p95") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% stats3a$redund.pred.names.p95]
    } else if (predictor.gr == "Rthresh_0p9_no_coords") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(stats3a$redund.pred.names.p90, params2$pred.names.sh$coords)]
    } else if (predictor.gr == "Rthresh_0p9_no_change") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(stats3a$redund.pred.names.p90, params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr)]
    } else if (predictor.gr == "Rthresh_0p9") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% stats3a$redund.pred.names.p90]
    }
    
    nr.vars <- length(predictors)
    
    if (params3a$mtry == 'sqrt_nr_var') {
      mtries <- floor(sqrt(nr.vars))
    } else if (params3a$mtry == 'nr_var_div_3') {
      mtries <- floor(nr.vars/3)
    }
  
    ## Random Forest
    RF.perf <- vector('double')
    for (targ in params3a$targ.names.lg) {
      
      if (!is.factor(X.trn$Ch_attr) | !is.factor(X.val$Ch_attr)) {
        stop("Ch_attr is not set as a factor")
      }
      
      ## Training
      if (params3a$parallel.RF) {
        set.seed(paramsGL$global.seed)
        nr.clusters <- min(params3a$ntree.MS, detectCores())
        tot.nrtrees <- params3a$ntree.MS
        nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)
        remaind <- tot.nrtrees%%(nr.clusters-1)
        cl <- makeCluster(nr.clusters)
        registerDoParallel(cl)
        rf.RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages='randomForest') %dopar% {
          randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=ntrees, mtry=mtries, nodesize=params3a$nodesize)
        }
        stopCluster(cl)
      } else { 
        set.seed(paramsGL$global.seed)
        rf.RF <- randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=params3a$ntree.MS, mtry=mtries, nodesize=params3a$nodesize)
      }
      
      ## Validation
      Y.val.predicted <- predict(rf.RF, X.val[,predictors], type="response", predict.all=F, nodes=F)
      if (!identical(rownames(X.val), rownames(Y.val))) {
        stop("TRN and VAL rownames do not match!!!")
      }
      cmd=sprintf('perf <- as.vector(regr_metrics(as.matrix(Y.val.predicted), Y.val$%s)[[params3a$metrics.MS]])', targ)   ## double [[]] to get only the value and not the name of list element
      eval(parse(text=cmd))
      RF.perf <- c(RF.perf, perf)
      
    }
    
    RF.entry.name <- sprintf('RF.%s.nrvars%s', predictor.gr, nr.vars)
    RF.entry <- as.list(c(RF.entry.name, RF.perf))
    list.RES_MS <- list(RES_MS.df, RF.entry)
    RES_MS.df <- rbindlist(list.RES_MS)
    
    ## YaImpute
    if (predictor.gr == params3a$best.predictor.group) {
      
      print('YaImpute')

      set.seed(paramsGL$global.seed) 
      
      ## - nodesize cannot be set
      ## - bootstrap=F to be used otherwise Error in yai(...: object 'xcvRefs' not found
      ## - rfMode="regression" to be used in all calls to yai(), otherwise it runs classification RF on discretized Ys (option "buildClasses") 
      ## - params3a$ntree.MS*length(params3a$targ.yaImp.names.lg) bc each Y uses a number of trees equal to ntree/nr.Ys
      yai.rf <- yai( x=X.trn[,predictors], y=Y.trn[,unlist(params3a$targ.yaImp.names.lg)], bootstrap=F, method="randomForest", rfMode="regression", k=1, ntree=params3a$ntree.MS*length(params3a$targ.yaImp.names.lg), mtry=mtries)
      
      newtargets.output <- newtargets(yai.rf, X.val[,predictors], k=1)
      ids.val.predicted <- newtargets.output$neiIdsTrgs
      ids.val.predicted <- cbind(rownames(ids.val.predicted), ids.val.predicted)
      colnames(ids.val.predicted) <- c("FCID", "Predicted_FCID")
      
      Y.trn.to.assign <- cbind(rownames(Y.trn), Y.trn)
      colnames(Y.trn.to.assign)[1] <- "FCID"
      
      ## join by FCID the two dataframes
      Y.val.predicted <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y = "FCID")
      rm(Y.trn.to.assign)
      
      ## sort both observed and predicted dataframes to have same row order
      Y.val.predicted <- arrange(Y.val.predicted, FCID)
      Y.val.observ <- cbind(rownames(Y.val), Y.val)
      colnames(Y.val.observ)[1] <- "FCID"
      Y.val.observ <- arrange(Y.val.observ, FCID)
      
      yaImp.perf <- vector('double')
      for (targ in params3a$targ.names.lg) { 
        cmd=sprintf('perf <- as.vector(regr_metrics(as.matrix(Y.val.predicted$%s), Y.val.observ$%s)[[params3a$metrics.MS]])', targ, targ)   ## double [[]] to get only the value and not the name of list element
        eval(parse(text=cmd))
        yaImp.perf <- c(yaImp.perf, perf)
      }
      rm(Y.val.observ)
      
      yaImp.entry.name <- sprintf('yaImp.%s.nrvars%s', params3a$best.predictor.group, nr.vars)
      yaImp.entry <- as.list(c(yaImp.entry.name, yaImp.perf))
      list.RES_MS <- list(RES_MS.df, yaImp.entry)
      RES_MS.df <- rbindlist(list.RES_MS)
      
    }
    
    temp.toc <- proc.time()-temp.tic[3]
    print(paste(predictor.gr, "elapsed time:", seconds_to_period(temp.toc[3])))
      
  }
  
  RES_MS.df <- t(RES_MS.df)
  write.csv(RES_MS.df, file = file.path(base_results_dir, sprintf("RES_MS_%s.csv", params3a$metrics.MS), sep = ''))

  stats3a$RES_MS.df <- RES_MS.df  ## add table to stats list to be saved
  
}

#### MODEL ASSESSMENT ---------------------------------------------------------------------------

if (params3a$run.MA) {

  ## selection only among plausible models
  if (params3a$best.predictor.group == "all") {
    predictors <- unlist(params2$pred.names.sh)
  } else if (params3a$best.predictor.group == "Rthresh_0p8") {
    predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% stats3a$redund.pred.names.p80]
  } else if (params3a$best.predictor.group == "Rthresh_0p95") {
    predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% stats3a$redund.pred.names.p95]
  } else if (params3a$best.predictor.group == "Rthresh_0p9") {
    predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% stats3a$redund.pred.names.p90]
  }
  
  nr.vars <- length(predictors)
  
  if (params3a$mtry == 'sqrt_nr_var') {
    mtries <- floor(sqrt(nr.vars))
  } else if (params3a$mtry == 'nr_var_div_3') {
    mtries <- floor(nr.vars/3)
  }
  
  ## Random Forest
  print('Model assessment')
  temp.tic <- proc.time() ## start clocking time for RF model assessment
  
  ## create model subdirectories to save models for later use (in mapping phase)
  models.subdir <- file.path(base_wkg_dir, "Models", fsep = .Platform$file.sep)
  if (! file.exists(models.subdir)){dir.create(models.subdir, showWarnings = F, recursive = T)}
  
  Y.val.predicted <- data.frame(matrix(nrow=nrow(X.val), ncol=length(params3a$targ.names.lg))) ## initialize matrix to store predictions
  rownames(Y.val.predicted) <- rownames(X.val)
  colnames(Y.val.predicted) <- params3a$targ.names.lg
  for (targ in params3a$targ.names.lg) {
    
    if (!is.factor(X.trn$Ch_attr) | !is.factor(X.val$Ch_attr)) {
      stop("Ch_attr is not set as a factor")
    }
    
    ## training
    if (params3a$parallel.RF) {
      set.seed(paramsGL$global.seed)
      nr.clusters <- min(params3a$ntree.MA, detectCores())
      tot.nrtrees <- params3a$ntree.MA
      nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)
      remaind <- tot.nrtrees%%(nr.clusters-1)
      cl <- makeCluster(nr.clusters)
      registerDoParallel(cl)
      rf.RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages='randomForest') %dopar% {
        randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=ntrees, mtry=mtries, nodesize=params3a$nodesize)
      }
      stopCluster(cl)
    } else {
      set.seed(paramsGL$global.seed)
      rf.RF <- randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=params3a$ntree.MA, mtry=mtries, nodesize=params3a$nodesize)
    }
    
    ## save RF model for each target
    cmd <- sprintf('RF.model.path <- file.path(models.subdir, "RF_%s.Rdata", fsep=.Platform$file.sep)', targ)   ## double [[]] to get only the value and not the name of list element
    eval(parse(text=cmd))
    save(rf.RF, file=RF.model.path)
    
    ## validation
    Y.val.predicted[,targ] <- predict(rf.RF, X.val[,predictors], type="response", predict.all=F, nodes=F)
    
  }
  
  ## save predictions to be loaded in case we want only to redo plots
  write.csv(Y.val.predicted, file = file.path(base_results_dir, sprintf("Y_val_predicted_matrix.csv"), sep = ''))

} else {
  Y.val.predicted <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_matrix.csv"), fsep = .Platform$file.sep))
}

#### STATS AND GRAPHS ---------------------------------------------------------------

## merge with initial datasets to get UTMzone, Ch_attr and ECOZONE columns for separate assessment wrt factors
Y.val.predicted <- rownames_2_FCID(Y.val.predicted)
Y.val.predicted <- merge(Y.val.predicted, X.trn.val[,c("FCID", "UTMzone", "Ch_attr")], by.x="FCID", by.y="FCID")
Y.val.predicted <- merge(Y.val.predicted, Y.trn.val[,c("FCID", "ECOZONE")], by.x="FCID", by.y="FCID")

## check order of samples in VAL set and rearrange dataframes with observed and predicted values to ensure that rows match
if (!identical(rownames(X.val), rownames(Y.val))) {
  stop("TRN and VAL rownames do not match!!!")
}
Y.val.predicted <- arrange(Y.val.predicted, FCID)
Y.val.observ <- rownames_2_FCID(Y.val)
Y.val.observ <- arrange(Y.val.observ, FCID)

## assessment
CAN.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.lg), ncol=length(params3a$metrics.MA.colnames) ) ) ## initialize matrix to store assessment stats (+2 bc range has 2 values)
colnames(CAN.stats.df) <- params3a$metrics.MA.colnames
rownames(CAN.stats.df) <- params3a$targ.names.lg 

ECO.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.lg), ncol=length(params3a$sampled.ecozones) ) ) 
colnames(ECO.stats.df) <- params3a$sampled.ecozones
rownames(ECO.stats.df) <- params3a$targ.names.lg

UTM.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.lg), ncol=length(paramsGL$zones) ) )
colnames(UTM.stats.df) <- paramsGL$zones
rownames(UTM.stats.df) <- params3a$targ.names.lg 

CNG.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.lg), ncol=length(params3a$Ch_attr.classes) ) )
colnames(CNG.stats.df) <- params3a$Ch_attr.classes
rownames(CNG.stats.df) <- params3a$targ.names.lg

titles <- c("mean height", "mean height", "mean height", "mean height", "mean height", "mean height", "mean height", "mean height", "mean height")
Assess.CAN.subdir <- file.path(base_figures_dir, "Assessment_CAN_level", fsep = .Platform$file.sep)
if (! file.exists(Assess.CAN.subdir)){dir.create(Assess.CAN.subdir, showWarnings = F, recursive = T)}

idx.row <- 1
for (targ in params3a$targ.names.lg) {
  
  cmd <- sprintf('temp.df <- data.frame(Y.val.predicted$%s, Y.val.observ$%s, Y.val.predicted[, c("ECOZONE", "UTMzone")], X.val$Ch_attr)', targ, targ)
  eval(parse(text=cmd))
  colnames(temp.df) <- c("predicted", "observed", "ECOZONE", "UTMzone", "Ch_attr")
  
  ## CAN-level stats
  CAN.stats.df[idx.row, ] <- as.vector(unlist(regr_metrics(temp.df$predicted, temp.df$observed)[params3a$metrics.MA]))
  
  ## CAN-level scatterplots
  str <- file.path(Assess.CAN.subdir, sprintf("ObsPred_scatter_%s.pdf", targ), sep='')
  pdf(str)
    scatt <- ggplot(temp.df, aes(x=observed, y=predicted)) +
      geom_point(shape=1) +    # Use hollow circles
      geom_abline(slope=1, intercept=0)  # Add 1:1 line line
     
    # XXXXXXXXXXXXXXXX 
#     - get lims
#     - set lims
    # XXXXXXXXXXXXXXXX
      coord_fixed(xlim = ) +
      ggtitle(sprintf("%s : R square=%.3f", titles[idx.row], CAN.stats.df[idx.row, "rsq"] )) +
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5))
  dev.off() 
  
  ## by Ecozone stats
  stats.by <- ddply(temp.df, "ECOZONE", summarize, regr_metrics(predicted, observed)[params3a$metrics.MS])   ## apply regr_metrics to the columns "predicted" and "observed" of dataframe temp.df split by "ECOZONE"...
  ECO.stats.df[idx.row, ] <- t(stats.by[match(as.character(params3a$sampled.ecozones), as.character(stats.by$ECOZONE)), ][,2])   ## ...then match the order of the stats.by dataframe to the list of Ecozones (in params3a) and finally take only 2nd column (with results) and transpose it to fill ECO.stats.df 
  
  ## by UTMZone stats
  zone.nrs <- data.frame("UTMzone"=substr(paramsGL$zones, 4, nchar(paramsGL$zones)))  ## to get rid of the "UTM" characters
  stats.by <- ddply(temp.df, "UTMzone", summarize, regr_metrics(predicted, observed)[params3a$metrics.MS])
  UTM.stats.df[idx.row, ] <- t(stats.by[match(as.character(zone.nrs$UTMzone), as.character(stats.by$UTMzone)), ][,2])
  
  ## by Ch_attr stats
  ############ TO DEL
#   if (any(table(temp.df$Ch_attr) <= 1)) {   ## to remove classes with 0 or 1 samples (error in lm() fit inside regr_metrics), happens only with subsampled data
#     ok.classes <- as.integer(levels(temp.df$Ch_attr)[table(temp.df$Ch_attr) > 1])
#     temp.df <- temp.df[temp.df$Ch_attr %in% ok.classes, ]
#     temp.df$Ch_attr <- factor(temp.df$Ch_attr)
#     warning("Some change types Ch_attr have only 0 or 1 samples in VAL set")
#   }
#   row.with.null <- t(stats.by[match(params3a$Ch_attr.labels, stats.by$Ch_attr), ][,2])
#   row.with.null[unlist(lapply(row.with.null, is.null))] <- NA
  ############ TO DEL
  
  stats.by <- ddply(temp.df, "Ch_attr", summarize, regr_metrics(predicted, observed)[params3a$metrics.MS])   
  CNG.stats.df[idx.row, ] <- t(stats.by[match(params3a$Ch_attr.labels, stats.by$Ch_attr), ][,2])

  
  idx.row <- idx.row+1
  
}


ECO.stats.df <- rbind(round(ECO.stats.df, 3), as.integer(TRN.samples.per.Ecozone), VAL.samples.per.Ecozone)
rownames(ECO.stats.df)[length(params3a$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
UTM.stats.df <- rbind(round(UTM.stats.df, 3), TRN.samples.per.UTMzone, VAL.samples.per.UTMzone)
rownames(UTM.stats.df)[length(params3a$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
CNG.stats.df <- rbind(round(CNG.stats.df, 3), TRN.samples.per.Ch_attr, VAL.samples.per.Ch_attr)
rownames(CNG.stats.df)[length(params3a$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")

stats3a$CAN.stats.df <- round(CAN.stats.df, 2)
stats3a$CAN.stats.df[, c("rsq", "nrmsd", "ac_uns", "ac_sys")] <- round(CAN.stats.df[, c("rsq", "nrmsd", "ac_uns", "ac_sys")], 3)
stats3a$ECO.stats.df <- ECO.stats.df
stats3a$UTM.stats.df <- UTM.stats.df
stats3a$CNG.stats.df <- CNG.stats.df

stats.file = file.path(base_results_dir, 'stats3a.Rdata', fsep = .Platform$file.sep) 
save(stats3a, file = stats.file)

write.csv(stats3a$CAN.stats.df, file = file.path(base_results_dir, "RES_MA_CAN_stats.csv", sep = ''))
write.csv(stats3a$ECO.stats.df, file = file.path(base_results_dir, sprintf("RES_MA_ECO_%s.csv", params3a$metrics.MS), sep = ''))
write.csv(stats3a$UTM.stats, file = file.path(base_results_dir, sprintf("RES_MA_UTM_%s.csv", params3a$metrics.MS), sep = ''))
write.csv(stats3a$CNG.stats.df, file = file.path(base_results_dir, sprintf("CNG_%s.csv", params3a$metrics.MS), sep = ''))

temp.toc <- proc.time()-temp.tic[3]
print( paste("Model assessment elapsed time:", seconds_to_period(temp.toc[3])) )


############# = TODEL
# yaImp.entry.name <- sprintf('yaImp.ntree%s.mtry%s', nrtrees, mtries)
# yaImp.entry <- as.list(c(yaImp.entry.name, yaImp.perf))
#
## variable importance plot
# yaiVarImp.500 <- yaiVarImp(yai.rf.500, nTop =0, plot=TRUE)
# str(yaiVarImp.500)
# par(las=1)
# par(mar=c(4,15,1,1))
#   boxplot(yaiVarImp.500, horizontal = TRUE, xlab = "Scaled Variable Importance")
# dev.off()
#
##save data and models as an RData file for future easy loading
## save rf and cca models
## save plot data both 6 response vars used in models and anncilary plot vars to be attached to maps
## save predcitor vars associated with plots
# save(yai.rf.500, X.trn, Y.trn.sub, Y.trn.attach, file = "yai_rf_gnn_sk.RData")
############# = TODEL

############# TOKEEPINCASE
## Not worth doing it in parallel for just 20'000 validation samples (to keep block of code in case we implement the mapping all over Canada in R)
#         ## Validation in parallel
#         sizeBlocks <- ceiling( nrow(X.val) / min(detectCores(), nrow(X.val)) )  ## to avoid problems when number of rows of the matrix is smaller than the nr of cores 
#         nblocks <- ceiling(nrow(X.val)/sizeBlocks)  ## to avoid problems when sizeBlocks is small (nblocks will allow to cover just a little more, or exactly, the size of data)
#         nr.clusters <- nblocks
#         cl <- makeCluster(nr.clusters)
#         registerDoParallel(cl)
#         Y.val.pred <- foreach (blck = 1:nblocks, .combine='rbind', .packages='randomForest') %do% {
#           # for (blck in 1:nblocks) {
#           indPart <- seq(from=(blck-1)*sizeBlocks+1, to=min(blck*sizeBlocks, nrow(X.val)), by=1)
#           Y.val.pred.block <- predict(rf, X.val[indPart,list.X.vars], type="response", predict.all=F, nodes=F)
#           as.data.frame(Y.val.pred.block)
#         }
#         stopCluster(cl)
#         cmd=sprintf('assess.rf.%s.val.ntree%s.mtry%s.nodesize%s <- regr_metrics(as.matrix(Y.val.pred), Y.val$%s)[list.metrics]', targ, nrtrees, mtries, nodesizes, targ)
#         eval(parse(text=cmd))
############# TOKEEPINCASE
       

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
print(paste("Prog3a, total elapsed time:",seconds_to_period(toc[3])))