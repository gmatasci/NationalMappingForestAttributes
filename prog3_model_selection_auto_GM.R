#### CODE INFOS -------------------------------------------------------------------

## Project Name: NationalImputationForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name:                            
## Objective: Descriptive statistics and plots, variable selection and Random Forest prediction/imputation

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - remove redefinition of base_wkg_dir
# - check final UTM zones to sample (now 17 with the removal of 11S)
# - params3a$parallel.RF.MS <- T 
# - params3a$parallel.RF.MA <- T (or F if we want RF importance)
# - params3a$subsetting <- F

# - reorder response variables to match order we want in paper
# - remove from training set aberrant samples wrt each response variable (these samples will not impact the predictions, producing values in a meaningful range --> RF, Yaimp both respect the range of training set)
# - ...or put in table 99.99% centile values
# -HAROLD: ask about ACuns and sys
# - IF PARALLEL RF IS USED SETTING SCALE=T WILL HAVE NO EFFECT, SO USE SEQUENTIAL RF FOR VARIABLE IMPORTANCE PLOT
# - run with final data: comparison between predicting response variables one at a time and not as a single multivariate Y: 
#   this means comparing yaimpute rf, rf, etc.
# - procedure to see sampled values based on "FCID" or "POLY250ID" from X or Y dataframes loaded here linked to field "POLY250ID" of <UTMzone>_pt_centerpt_training_validation.shp (saved in /BAP_Imputation_working/wkg/<UTMzone>/)
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
# -V check plots with new definition of year
# -V check content and order/sorting with arrange of different Y.val.predicted: RF vs. Yai
# -V check all lines where removal of Bands has an impact
# -V check yai object to see if 6 RF
# -v add variable importance plot for both (boxplot over the 10 or 6 values?) -- will not add Yai one bc it is not clear what is it that they plot (lack of info in docs)
# -V remove bias_pct from plots, add RMSE% (also in table, and remove nRMSE)
# -V fix empty csvs issue in the end
# -V produce barplots for ECO, UTM, CNG assessments

#### READS/WRITES ------------------------------------------------------------

## READS:
# - "lidar_metrics_mean_training_validation.csv": (from prog1) csv table with all observed LiDAR and forest attributes (Y) for selected TRN and VAL samples (3x3 polygons with average plot values or just single plot value) for all UTM zones
# - "poly_training_validation_exvars_extract.csv":  (from prog2) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

## WRITES:
# - "yai_rf_gnn_sk.RData": R data file with yaimpute models and training set to use later for mapping (prediction/imputation) on the grid

#### INIT --------------------------------------------------------------------

print('Prog3a: descriptive stats/plots, model selection and model assessment') 

rm(list=ls())

param_file <- "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

source("D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/Functions_NatImp.R")

#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params3a <- list()

## Actual parameters to be used
params3a$methods <- list("RF", "YAI")
params3a$run.checks <- F   ## where to run checks for duplicate FCIDs, run only once, as of 06/06/2016 all is OK
params3a$run.descr.stats <- F  ## whether to run or not descriptive stats block
params3a$run.MS <- F
params3a$run.MA <- F
params3a$parallel.RF.MS <- T
params3a$parallel.RF.MA <- F
params3a$subsetting <- F
params3a$ntree.MS <- 32   ## ideally 500 but then YaImp (not parallel) takes ages
params3a$ntree.MA <- 100

## Parameters for testing
# params3a$methods <- list("RF", "YAI")
# params3a$run.checks <- F
# params3a$run.descr.stats <- T
# params3a$run.MS <- T
# params3a$run.MA <- T
# params3a$parallel.RF.MS <- T
# params3a$parallel.RF.MA <- F
# params3a$subsetting <- T
# params3a$ntree.MS <- 10
# params3a$ntree.MA <- 20

params3a$mtry <- 'sqrt_nr_var'
params3a$nodesize <- 5
params3a$plot.importance <- T

# params3a$predictor.groups <- list("all", "SKstudy", "Rthresh_0p8", "Rthresh_0p95", "Rthresh_0p9_no_coords", "Rthresh_0p9_no_change", "Rthresh_0p9")
params3a$predictor.groups <- list("all", "Rthresh_0p9")
params3a$best.predictor.group <- "Rthresh_0p9"

params3a$predictors.to.rem.prior.know <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7')  ## do not consider all of the raw spectral bands (only removing VIS bands here results in the same set as b4, b5 and b7 are then filtered out based on |R|)

params3a$metrics.MS <- "rsq"
params3a$metrics.MA <- c("xmean", "xrangemin", "xrangemax", "xstdev", "ymean", "yrangemin", "yrangemax", "ystdev", "rsq", "rmsd", "rmsdpct", "ac_uns", "ac_sys", "bias")
params3a$metrics.MA.colnames <- c("o_mean", "o_rangemin", "o_rangemax", "o_stdev", "p_mean", "p_rangemin", "p_rangemax", "p_stdev", "rsq", "rmsd", "rmsdpct", "ac_uns", "ac_sys", "bias")

## final target variables to be predicted by both RF and YaImp
params3a$targ.names.lg <- list("elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean", "loreys_height", "basal_area", "gross_stem_volume", "total_biomass")
params3a$targ.names.sh <- list("el_m","el_std","el_p95","el_cv","pct_1r_ab2","pct_1r_abm", "loreys_h", "basal_a", "stem_v", "tot_biom")

## target variables given as input to YaImp to find best neighbor (other 4 are attached once FCIDs are known) 
params3a$targ.yaImp.names.lg <- list("elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean")
params3a$targ.yaImp.names.sh <- list("el_m","el_std","el_p95","el_cv","pct_1r_ab2","pct_1r_abm")

## two key target variables used in descriptive plots to show X to Y relationship
params3a$targ.key.names.lg <- list("percentage_first_returns_above_2m", "total_biomass")
params3a$targ.key.names.sh <- list("pct_1r_ab2", "tot_biom")

## final target variable names and units to be used in plots
params3a$targ.names.plots <- list("elev_mean", "elev_sd", "elev_p95", "elev_cv", "cover_2m", "cover_mean", "loreys_height", "basal_area", "stem_volume", "ag_biomass")
params3a$targ.units.plots <- list("m", "m", "m", "", "%", "%", "m", "m^2/ha", "m^3/ha", "kg/ha")

## retained 5 classes: 0="No change", 1="Fire", 2="Harvesting", 3="Non-stand replacing", 4="Infrastructure"
params3a$Ch_attr.classes <- c("No change", "Fire", "Harvesting", "Non-stand replacing", "Infrastructure")
params3a$Ch_attr.labels <- c(0, 1, 2, 3)  ## no "Infrastructure" bc there are too few points

params3a$groups.to.plot <- list('Bands', 'TCcomps_VI', 'ChangeAttribution', 'Years_Topo_Coords')  ## to plot relation with key response vars separately (and in a visually pleasing way)

params3a$nr.pts.plot <- 2000 ## to make pairs plots less dense in points and to subset the dataset for tests in the developement phase

params3a$corr.high <- 0.9  ## corr threshold to highlight high correlation among predictors
params3a$corr.low <- 1  ## corr threshold to highlight low correlation between predictors and response vars

params3a$nr.of.bins <- 12  ## nr of bins for UTM zone histograms
params3a$hist.centile <- 0.999  ## percentile threshold to subset data for histograms (avoid extreme values setting a too high xlim)

params3a$all.ecozones <- list("Arctic Cordillera", "Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", 
                          "Hudson Plains", "Mixedwood Plains", "Montane Cordillera", "Northern Arctic", "Pacific Maritime", "Semiarid Prairies", 
                          "Southern Arctic", "Subhumid Prairies", "Taiga Cordillera", "Taiga Plains", "Taiga Shield East", "Taiga Shield West")

params3a$sampled.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", "Hudson Plains", 
                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")  ## Semiarid Prairies and Subhumid Prairies (ECOZONE_ID = 10 in shp) result in 0 and 3 plots, respectively, so we remove them

param_file_prog3a = file.path(base_wkg_dir, 'AllUTMzones_params3a.Rdata', fsep = .Platform$file.sep) 
save(params3a, file = param_file_prog3a)

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("latex2exp",
                      "xtable",
                      "rgeos",
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
idx.raster.zero <- idx.Landsat.zero | idx.Elev.zero  ## set as TRUE if either of the 2 conditions is TRUE
Y.trn.val <- Y.trn.val.raw[!idx.raster.zero, ]  ## keep rows that are not (!) to remove
X.trn.val <- X.trn.val.raw[!idx.raster.zero, ]
rm(Y.trn.val.raw, X.trn.val.raw)  ## clear space after filtering

## delete samples from Subhumid Prairies (just 3 plots after filtering)
idx.Subhumid.Prairies <- Y.trn.val$ECOZONE == "Subhumid Prairies"
Y.trn.val <- Y.trn.val[!idx.Subhumid.Prairies, ]
X.trn.val <- X.trn.val[!idx.Subhumid.Prairies, ]
Y.trn.val$ECOZONE <- factor(Y.trn.val$ECOZONE)  ## to drop the level we just removed

## save stats about removed samples and nr of sample from each Ecozone, UTM zone and Ch_attr
stats3a <- list()
stats3a$NA.samples.removed <- sum(idx.raster.zero==T) + sum(idx.Subhumid.Prairies==T)
TRN.samples.per.Ecozone <- table(Y.trn.val[Y.trn.val$TV=="TRAINING",]$ECOZONE)  ## to be added later on to the final stats dataframes
VAL.samples.per.Ecozone <- table(Y.trn.val[Y.trn.val$TV=="VALIDATION",]$ECOZONE)
TRN.samples.per.UTMzone <- table(X.trn.val[X.trn.val$TV=="TRAINING",]$UTMzone)
VAL.samples.per.UTMzone <- table(X.trn.val[X.trn.val$TV=="VALIDATION",]$UTMzone)
TRN.samples.per.Ch_attr <- table(X.trn.val[X.trn.val$TV=="TRAINING",]$Ch_attr)
VAL.samples.per.Ch_attr <- table(X.trn.val[X.trn.val$TV=="VALIDATION",]$Ch_attr)

#### CORRELATION MATRIX ----------------------------------------------

## compute the following outside of if block bc we want to know redundant predictors in all cases
prior.know.cont.idx <- !unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr) ## remove unwanted predictors by prior knowledge and not-continuous predictors
X.trn.val.corr <- X.trn.val[, unlist(params2$pred.names.sh)[prior.know.cont.idx]]
X.corr.matrix <- cor(X.trn.val.corr)  ## compute corr only among predictors (top-left sub-block)

## print correlation matrix as a Latex table
print( xtable(X.corr.matrix, digits=2), file = file.path(base_results_dir, "X_corr_matrix.tex", sep = '') )

## based on the dataset after removal of unwanted predictors by prior knowledge, find subset of redundant variables automatically:
## in a pair having r>cutoff remove the variable with highest average absolute correlation with rest of available variables (exact=T)
stats3a$redund.pred.names.p80 <- findCorrelation(X.corr.matrix, cutoff=0.80, exact=T, verbose=F, names=T)
stats3a$redund.pred.names.p90 <- findCorrelation(X.corr.matrix, cutoff=0.90, exact=T, verbose=F, names=T)
stats3a$redund.pred.names.p95 <- findCorrelation(X.corr.matrix, cutoff=0.95, exact=T, verbose=F, names=T)

if (params3a$run.descr.stats) {
  
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plot2 <- ggplot(df.to.plot, aes(x=Ch_attr, y=tot_biom, fill=Ch_attr)) + 
        geom_boxplot(notch=F, outlier.shape = NA) + 
        scale_x_discrete(labels=params3a$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        coord_cartesian(ylim = c(0, 250000))
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
  cat.vars <- c("YrsSince_GrCh", "Ch_attr")  ## categorical vars to remove
  cont.var.idx <- !(colnames(X.trn.val) %in% cat.vars)  ## get indixes of all predictors except "Ch_attr" and "YrsSince_GrCh" (dealt separately as a barplot)
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
      if (var == "YrsSince_GrCh") {
        yr.range <- 0:25  ## full range to be barplotted
        yr.padded <- c( data.to.plot[zone.idx, 6], yr.range[!yr.range %in% unique(data.to.plot[zone.idx, 6])] )  ## add to the year data the years that are missing to have a complete set 
        barplot(table(yr.padded[yr.padded!=50]), main=sprintf("%s", zone), ylim=c(0, 300))  ## remove no-change data, years never appearing will have a value of 1 instead of 0 in the bar plot 
      } else {
        barplot(table(data.to.plot[zone.idx, 6]), main=sprintf("%s", zone), ylim=c(0, 1200))
      }
    }
    title(var, outer=T)
    dev.off()
  }

}

#### DATA PREPARATION -------------------------------------------------

## remove unneeded vars and predictors that we discard by a priori knowledge in X.trn.val data
stats3a$initial.predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% params3a$predictors.to.rem.prior.know]  ## save list of predictors at the start of the analysis (after removal of a priori exclusions) 
X.trn.val.data <- X.trn.val[, c('FCID', stats3a$initial.predictors, 'TV')]

## extract core str response variables using in random forest analysis (leave out ECOZONE, to be retrieved later on through the FCID)
## imputation map only contains single band of plot FCIDs bc extracted key variables that get attached to map after imputation
Y.trn.val.data <- Y.trn.val[, c("FCID", unlist(params3a$targ.names.lg), "TV")]

X.trn.val.data <- FCID_2_rownames(X.trn.val.data)
Y.trn.val.data <- FCID_2_rownames(Y.trn.val.data)

if (params3a$subsetting) {
  set.seed(paramsGL$global.seed)
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
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% params3a$predictors.to.rem.prior.know]  ## always remove the a priori excluded variables
    } else if (predictor.gr == "SKstudy") {
      predictors <- c("TCB", "TCG", "TCW", "TCA", "TCD", "YrsSince_GrCh", "Ch_attr", "Elev", "Slope", "TWI", "TSRI", "Long", "Lat")  ## 13 predictors
    } else if (predictor.gr == "Rthresh_0p8") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, stats3a$redund.pred.names.p80)]
    } else if (predictor.gr == "Rthresh_0p95") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, stats3a$redund.pred.names.p95)]
    } else if (predictor.gr == "Rthresh_0p9_no_coords") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, stats3a$redund.pred.names.p90, params2$pred.names.sh$coords)]
    } else if (predictor.gr == "Rthresh_0p9_no_change") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, stats3a$redund.pred.names.p90, params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr)]
    } else if (predictor.gr == "Rthresh_0p9") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, stats3a$redund.pred.names.p90)]
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
      if (params3a$parallel.RF.MS) {
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
      cmd=sprintf('perf <- as.vector(regr_metrics(Y.val$%s, as.matrix(Y.val.predicted))[[params3a$metrics.MS]])', targ)   ## double [[]] to get only the value and not the name of list element
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
      ids.val.predicted <- rownames_2_FCID(ids.val.predicted)
      colnames(ids.val.predicted)[2] <- "Predicted_FCID"
      
      Y.trn.to.assign <- rownames_2_FCID(Y.trn)
      
      ## join by FCID the two dataframes
      Y.val.predicted <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y = "FCID")
      rm(Y.trn.to.assign)
      
      ## sort both observed and predicted dataframes to have same row order
      Y.val.predicted <- arrange(Y.val.predicted, FCID)
      Y.val.observ <- rownames_2_FCID(Y.val)
      Y.val.observ <- arrange(Y.val.observ, FCID)
      
      yaImp.perf <- vector('double')
      for (targ in params3a$targ.names.lg) { 
        cmd=sprintf('perf <- as.vector(regr_metrics(Y.val.observ$%s, as.matrix(Y.val.predicted$%s))[[params3a$metrics.MS]])', targ, targ)   ## double [[]] to get only the value and not the name of list element
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

## create model subdirectories to save models for later use (in mapping phase)
models.subdir <- file.path(base_wkg_dir, "Models", fsep = .Platform$file.sep)
if (! file.exists(models.subdir)){dir.create(models.subdir, showWarnings = F, recursive = T)}

## selection only among plausible models
if (params3a$best.predictor.group == "all") {
  predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% params3a$predictors.to.rem.prior.know]
} else if (params3a$best.predictor.group == "Rthresh_0p9") {
  predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3a$predictors.to.rem.prior.know, stats3a$redund.pred.names.p90)]
}

nr.vars <- length(predictors)

stats3a$final.predictors <- predictors

if (params3a$run.MA) {
  
  print('Model assessment')
  
  if (params3a$mtry == 'sqrt_nr_var') {
    mtries <- floor(sqrt(nr.vars))
  } else if (params3a$mtry == 'nr_var_div_3') {
    mtries <- floor(nr.vars/3)
  }

#### RANDOM FOREST ---------------------------------------------------------------------
  
  print('Random Forest')
  temp.tic <- proc.time()
  
  Y.val.predicted.RF <- data.frame(matrix(nrow=nrow(X.val), ncol=length(params3a$targ.names.lg))) ## initialize matrix to store predictions
  rownames(Y.val.predicted.RF) <- rownames(X.val)
  colnames(Y.val.predicted.RF) <- params3a$targ.names.lg
  Y.val.predicted.stddev.RF <- Y.val.predicted.RF
  for (targ in params3a$targ.names.lg) {
    
    if (!is.factor(X.trn$Ch_attr) | !is.factor(X.val$Ch_attr)) {
      stop("Ch_attr is not set as a factor")
    }
    
    ## Training
    if (params3a$parallel.RF.MA) {
      set.seed(paramsGL$global.seed)
      nr.clusters <- min(params3a$ntree.MA, detectCores())
      tot.nrtrees <- params3a$ntree.MA
      nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)
      remaind <- tot.nrtrees%%(nr.clusters-1)
      cl <- makeCluster(nr.clusters)
      registerDoParallel(cl)
      rf.RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages='randomForest') %dopar% {
        randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=ntrees, mtry=mtries, nodesize=params3a$nodesize, importance=params3a$plot.importance)
      }
      stopCluster(cl)
    } else {
      set.seed(paramsGL$global.seed)
      rf.RF <- randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=params3a$ntree.MA, mtry=mtries, nodesize=params3a$nodesize, importance=params3a$plot.importance)
    }
    
    ## save RF model for each target
    cmd <- sprintf('RF.model.path <- file.path(models.subdir, "RF_%s.Rdata", fsep=.Platform$file.sep)', targ)
    eval(parse(text=cmd))
    save(rf.RF, file=RF.model.path)
    
    ## validation
    prediction.res <- predict(rf.RF, X.val[,predictors], type="response", predict.all=T, nodes=F)
    Y.val.predicted.RF[,targ] <- prediction.res$aggregate
    Y.val.predicted.stddev.RF[,targ] <- apply(prediction.res$individual, 1, sd)  ## std dev of each row
    
  }
  
  ## arrange and save predictions and uncertainties to be loaded in case we want only to redo plots
  Y.val.predicted.RF <- rownames_2_FCID(Y.val.predicted.RF)  ## to make the column FCID available for the merge() function in the assessment part below
  Y.val.predicted.RF <- arrange(Y.val.predicted.RF, FCID)
  Y.val.predicted.stddev.RF <- rownames_2_FCID(Y.val.predicted.stddev.RF)  ## to make the column FCID available for the merge() function in the assessment part below
  Y.val.predicted.stddev.RF <- arrange(Y.val.predicted.stddev.RF, FCID)
  write.csv(Y.val.predicted.RF, file = file.path(base_results_dir, sprintf("Y_val_predicted_matrix_RF.csv"), sep = ''))
  write.csv(Y.val.predicted.stddev.RF, file = file.path(base_results_dir, sprintf("Y_val_predicted_stddev_matrix_RF.csv"), sep = ''))
  
  temp.toc <- proc.time()-temp.tic[3]
  print(paste("Random Forest elapsed time:", seconds_to_period(temp.toc[3])))

#### YAIMPUTE -------------------------------------------------------------------
  
  print('YaImpute')
  temp.tic <- proc.time()
  
  set.seed(paramsGL$global.seed) 
  ## - nodesize cannot be set
  ## - bootstrap=F to be used otherwise Error in yai(...: object 'xcvRefs' not found
  ## - rfMode="regression" to be used in all calls to yai(), otherwise it runs classification RF on discretized Ys (option "buildClasses") 
  ## - params3a$ntree.MS*length(params3a$targ.yaImp.names.lg) bc each Y uses a number of trees equal to ntree/nr.Ys
  yai.rf <- yai( x=X.trn[,predictors], y=Y.trn[,unlist(params3a$targ.yaImp.names.lg)], bootstrap=F, method="randomForest", rfMode="regression", k=1, ntree=params3a$ntree.MA*length(params3a$targ.yaImp.names.lg), mtry=mtries)
  
  ## NOT CLEAR HOW IT IS COMPUTED
  ## variable importance plot
  # Yai.var.imp <- yaiVarImp(yai.rf, nTop =0, plot=F)
  # par(las=1)
  # par(mar=c(4,15,1,1))
  #   boxplot(Yai.var.imp, horizontal = TRUE, xlab = "Scaled Variable Importance")
  # dev.off()
  ## NOT CLEAR HOW IT IS COMPUTED
  
  ## impute IDs
  newtargets.output <- newtargets(yai.rf, X.val[,predictors], k=1)
  ids.val.predicted <- newtargets.output$neiIdsTrgs
  ids.val.predicted <- rownames_2_FCID(ids.val.predicted)
  colnames(ids.val.predicted)[2] <- "Predicted_FCID"
  
  Y.trn.to.assign <- rownames_2_FCID(Y.trn)
  
  ## join by FCID the two dataframes to get imputed values
  Y.val.predicted.Yai <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
  Y.val.predicted.Yai$Predicted_FCID <- NULL
  
  ## save model and Y trn to impute values when mapping
  Yai.model.path <- file.path(models.subdir, "YAI.Rdata", fsep=.Platform$file.sep)   
  save(yai.rf, file=Yai.model.path)
  Y.trn.path <- file.path(models.subdir, "Ytrn.Rdata", fsep=.Platform$file.sep)
  save(Y.trn.to.assign, file=Y.trn.path)
  rm(Y.trn.to.assign)
  
  ## arrange and save predictions to be loaded in case we want only to redo plots
  Y.val.predicted.Yai <- arrange(Y.val.predicted.Yai, FCID)
  write.csv(Y.val.predicted.Yai, file = file.path(base_results_dir, sprintf("Y_val_predicted_matrix_YAI.csv"), sep = ''))
  
  temp.toc <- proc.time()-temp.tic[3]
  print(paste("YaImpute elapsed time:", seconds_to_period(temp.toc[3])))

} else {  ## if we did not run the model assessment, read the csv with the predictions that we saved beforehand
  Y.val.predicted.RF <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_matrix_RF.csv"), fsep = .Platform$file.sep), row.names=1)  ## row.names=1 to have same format as when these dataframes are produced in the predictions phase
  Y.val.predicted.stddev.RF <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_stddev_matrix_RF.csv"), fsep = .Platform$file.sep), row.names=1)
  Y.val.predicted.Yai <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_matrix_YAI.csv"), fsep = .Platform$file.sep), row.names=1)
}

#### STATS AND GRAPHS ---------------------------------------------------------------

## subdirectory to save model assessment plots
Assess.CAN.subdir <- file.path(base_figures_dir, "Assessment_CAN_level", fsep = .Platform$file.sep)
if (! file.exists(Assess.CAN.subdir)){dir.create(Assess.CAN.subdir, showWarnings = F, recursive = T)}

for (method in params3a$methods) {

  if (method == "RF") {
    Y.val.predicted <- Y.val.predicted.RF
  } else if (method == "YAI") {
    Y.val.predicted <- Y.val.predicted.Yai
  }
  
  ## merge with initial datasets to get UTMzone, Ch_attr and ECOZONE columns for separate assessment wrt factors
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
  CAN.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.plots), ncol=length(params3a$metrics.MA.colnames) ) ) ## initialize matrix to store assessment stats (+2 bc range has 2 values)
  colnames(CAN.stats.df) <- params3a$metrics.MA.colnames
  rownames(CAN.stats.df) <- params3a$targ.names.plots 
  
  ECO.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.plots), ncol=length(params3a$sampled.ecozones) ) ) 
  colnames(ECO.stats.df) <- params3a$sampled.ecozones
  rownames(ECO.stats.df) <- params3a$targ.names.plots
  
  UTM.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.plots), ncol=length(paramsGL$zones) ) )
  colnames(UTM.stats.df) <- paramsGL$zones
  rownames(UTM.stats.df) <- params3a$targ.names.plots 
  
  CNG.stats.df <- data.frame( matrix( nrow=length(params3a$targ.names.plots), ncol=length(params3a$Ch_attr.classes)-1 ) )  ## we do not consider the class "Infrastructure" in the separate assessment by Ch_attr
  colnames(CNG.stats.df) <- params3a$Ch_attr.classes[-length(params3a$Ch_attr.classes)]
  rownames(CNG.stats.df) <- params3a$targ.names.plots
  
  plot.xlabel <- "observed"
  plot.ylabel <- "predicted"
  
  ## initialize empty dataframe to host the complete results table for shp: 10 targ * 3 values (observeed, predicted, residual) = 30 columns
  results.shp.df <- data.frame( matrix( nrow=nrow(Y.val.observ), ncol=3*length(params3a$targ.names.lg) ) ) 
  
  ## initialize variable importance dataframe
  var.imp.df <- data.frame( matrix( nrow=length(params3a$targ.names.lg), ncol=nr.vars ) ) ## initialize matrix to store variable importance values
  colnames(var.imp.df) <- predictors
  rownames(var.imp.df) <- params3a$targ.names.lg
  
  idx.targ <- 1
  for (targ in params3a$targ.names.lg) {
    
    ## load RF model for each target
    cmd <- sprintf('load(file.path(models.subdir, "RF_%s.Rdata", fsep=.Platform$file.sep))', targ) 
    eval(parse(text=cmd))

    ## fill variable importance dataframe
    ## type=1: mean decrease in accuracy (mean increase in MSE for regression). For each tree, the prediction error on the out-of-bag portion of the data is
    ## recorded. Then the same is done after permuting each predictor variable. The difference between the two are then averaged over all trees, and normalized by the standard deviation of the differences.
    ## IF PARALLEL RF IS USED SETTING SCALE=T WILL HAVE NO EFFECT, SO USE SEQUENTIAL RF FOR VARIABLE IMPORTANCE PLOT
    var.imp.df[targ, ] <- t(importance(rf.RF, type=1, scale=T))  
    
    ## build temporary df to ease by category stats with ddply()
    cmd <- sprintf('temp.df <- data.frame(Y.val.observ$%s, Y.val.predicted$%s, Y.val.predicted[, c("ECOZONE", "UTMzone")], X.val[, c("Ch_attr", "Long", "Lat")])', targ, targ)
    eval(parse(text=cmd))
    colnames(temp.df) <- c("observed", "predicted", "ECOZONE", "UTMzone", "Ch_attr", "Long", "Lat")
    
    ## store in results.shp.df all predictions and residuals to be written in a shp
    results.shp.df[, seq((idx.targ-1)*3+1, (idx.targ-1)*3+3)] <- data.frame(temp.df$observed, temp.df$predicted, temp.df$predicted - temp.df$observed)
    cmd <- sprintf('colnames(results.shp.df)[seq((idx.targ-1)*3+1, (idx.targ-1)*3+3)] <- c(\"O_%s\", \"P_%s\", \"R_%s\")', params3a$targ.names.sh[idx.targ],  params3a$targ.names.sh[idx.targ],  params3a$targ.names.sh[idx.targ])
    eval(parse(text=cmd))

    ## CAN-level stats
    CAN.stats.compl <- unlist(regr_metrics(temp.df$observed, temp.df$predicted)[c(params3a$metrics.MA, "b_yvsx", "a_yvsx")])
    CAN.stats.df[idx.targ, ] <- as.vector(CAN.stats.compl[params3a$metrics.MA])
      
    ## CAN-level scatterplots
    fig.name.str <- file.path(Assess.CAN.subdir, sprintf("ObsPred_scatter_%s_%s.pdf", method, targ), sep='')
    title.str <- sprintf("%s [%s]: ", params3a$targ.names.plots[idx.targ], params3a$targ.units.plots[idx.targ])
    box.str <- sprintf("R^2=%.3f, RMSE%%=%.1f, bias=%.2f", CAN.stats.df[idx.targ, "rsq"], CAN.stats.df[idx.targ, "rmsdpct"], CAN.stats.df[idx.targ, "bias"])
    scatt.range <- c(min(temp.df$observed, temp.df$predicted), quantile(c(temp.df$observed, temp.df$predicted), 0.9999, names=F))
    GMFR.slope <- CAN.stats.compl["b_yvsx"]
    GMFR.interc <- CAN.stats.compl["a_yvsx"]
    lm.mod <- lm(predicted~observed, data=temp.df)
    lm.slope <- lm.mod$coefficients[2]
    lm.interc <- lm.mod$coefficients[1]
    ## check if corr. coeff. R is equal to the slope of the regression based on standardized predicted and observed values: pred = a + b obs
    temp.df.stand <- as.data.frame(scale(temp.df[, c("predicted", "observed")], center=T, scale=T))
    lm.mod.stand <- lm(predicted~observed, data=temp.df.stand)
    if ( as.vector(lm.mod.stand$coefficients[2]) - as.vector(sqrt(CAN.stats.df[idx.targ, "rsq"])) > 1e-10 ) {
      warning(sprintf("Predicted VS observed R and standardized regression slope do not match for %s on %s", method, targ))
    }
    pdf(fig.name.str)
      theme_set(theme_gray(base_size = 18))
      scatt <- plot_colorByDensity(temp.df$observed, temp.df$predicted, xlim=scatt.range, ylim=scatt.range, xlab=plot.xlabel, ylab=plot.ylabel, main=paste(title.str, '\n', box.str)) +
                  geom_abline(intercept = GMFR.interc, slope = GMFR.slope, linetype="dashed", size=1) +   ## Add GMFR line
                  # geom_abline(intercept = lm.interc, slope = lm.slope, linetype="dashed", size=1, color='red') +   ## Add lm line
                  geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line
    print(scatt)
    dev.off()

    ## by Ecozone stats
    stats.by <- ddply(temp.df, "ECOZONE", summarize, regr_metrics(observed, predicted)[params3a$metrics.MS])   ## apply regr_metrics to the columns "predicted" and "observed" of dataframe temp.df split by "ECOZONE"...
    ECO.stats.df[idx.targ, ] <- t(stats.by[match(as.character(params3a$sampled.ecozones), as.character(stats.by$ECOZONE)), ][,2])   ## ...then match the order of the stats.by dataframe to the list of Ecozones (in params3a) and finally take only 2nd column (with results) and transpose it to fill ECO.stats.df 
    
    ## by UTMZone stats
    zone.nrs <- data.frame("UTMzone"=substr(paramsGL$zones, 4, nchar(paramsGL$zones)))  ## to get rid of the "UTM" characters
    stats.by <- ddply(temp.df, "UTMzone", summarize, regr_metrics(observed, predicted)[params3a$metrics.MS])
    UTM.stats.df[idx.targ, ] <- t(stats.by[match(as.character(zone.nrs$UTMzone), as.character(stats.by$UTMzone)), ][,2])
    
    ## by Ch_attr stats
    temp.df.Ch <- temp.df[temp.df$Ch_attr!=4, ]
    temp.df.Ch$Ch_attr <- factor(temp.df.Ch$Ch_attr)  ## to drop the level 4
    stats.by <- ddply(temp.df.Ch, "Ch_attr", summarize, regr_metrics(observed, predicted)[params3a$metrics.MS])   
    CNG.stats.df[idx.targ, ] <- t(stats.by[match(params3a$Ch_attr.labels, stats.by$Ch_attr), ][,2])
    
    idx.targ <- idx.targ+1
    
  }
  
  ## variable importance plot
  if (params3a$plot.importance) {
    
    rownames(var.imp.df) <- params3a$targ.names.plots
    
    ## overall boxplot over the 10 targets
    fig.name.str <- file.path(Assess.CAN.subdir, sprintf("RF_VariableImportance_IncMSE_overall.pdf"), sep='')
    pdf(fig.name.str)
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 4, 0, 0)) 
      theme_set(theme_gray(base_size = 18))
      medians <- sapply(var.imp.df, median)
      idx.medians <- order(medians, decreasing=T)
      boxpl <- boxplot(var.imp.df[, idx.medians], horizontal = TRUE, xlab = "Scaled Variable Importance (% incr. MSE)", las=2)
    print(boxpl)
    dev.off()
    
    ## heatmap by target
    fig.name.str <- file.path(Assess.CAN.subdir, sprintf("RF_VariableImportance_IncMSE_byTarg.pdf"), sep='')
    pdf(fig.name.str)
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 4, 0, 0)) 
      theme_set(theme_gray(base_size = 18))
      my_palette <- colorRampPalette(c("white", "dark green"))(n = 20)
      var.imp.df.to.plot <- round(as.matrix(var.imp.df[, idx.medians]), 1)  ## sort by descending median value
      heatmap.2(var.imp.df.to.plot,
                cellnote = var.imp.df.to.plot,  # same data set for cell labels
                notecol="black",      # change font color of cell labels to black
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(9,9),
                col=my_palette,       # use on color palette defined earlier
                key=TRUE, 
                key.title="", 
                key.xlab="Scaled Variable Importance (% incr. MSE)",
                lmat=rbind(c(2,3), c(0,1), c(0,4)),
                lwid=c(1.5, 4),
                lhei=c(1.5, 4, 1.2),
                dendrogram="none",     # only draw a row dendrogram
                Rowv=F,
                Colv=F
                )   
    print(boxpl)
    dev.off()
    
  }
  
  ## write shp with individual validation plot prediction results
  results.shp.df <- data.frame(results.shp.df, X.val$Long, X.val$Lat)
  coordinates(results.shp.df) <- ~X.val.Long+X.val.Lat
  proj4string(results.shp.df) <- CRS("+proj=longlat +datum=NAD83")
  writeOGR(results.shp.df, base_results_dir, sprintf('val_predictions_%s', method), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## just to check validity of other layers, never used later on
  
  ## correlation and corr difference matrices for observed/predicted values on validation set
  if (method == 'RF') {
    obs.data.matrix <- results.shp.df@data[ , seq(1,ncol(results.shp.df),3)]
    colnames(obs.data.matrix) <- params3a$targ.names.plots
    obs.corr.matrix <- cor(obs.data.matrix)
    write.csv(obs.corr.matrix, file = file.path(base_results_dir, "Val_obs_corr_matrix.csv", sep = ''))
    print( xtable(obs.corr.matrix, digits=rep(2, length(ncol(obs.corr.matrix)))), file = file.path(base_results_dir, "Val_obs_corr_matrix.tex", sep = '') )
  }
  pred.data.matrix <- results.shp.df@data[ , seq(2,ncol(results.shp.df),3)]
  colnames(pred.data.matrix) <- params3a$targ.names.plots
  pred.corr.matrix <- cor(pred.data.matrix)
  write.csv(pred.corr.matrix, file = file.path(base_results_dir, sprintf("Val_pred_%s_corr_matrix.csv", method), sep = ''))
  cmd <- sprintf( "print( xtable(pred.corr.matrix, digits=rep(2, length(ncol(pred.corr.matrix)))), file = file.path(base_results_dir, \"Val_pred_%s_corr_matrix.tex\", sep = '') )", method)
  eval(parse(text=cmd))
  diff.corr.matrix <- pred.corr.matrix-obs.corr.matrix
  write.csv(diff.corr.matrix, file = file.path(base_results_dir, sprintf("Val_diff_%s_corr_matrix.csv", method), sep = ''))
  cmd <- sprintf( "print( xtable(diff.corr.matrix, digits=rep(2, length(ncol(diff.corr.matrix)))), file = file.path(base_results_dir, \"Val_diff_%s_corr_matrix.tex\", sep = '') )", method)
  eval(parse(text=cmd))
  
  
  ## TOO COMPLICATED
  # ## Maps of residuals
  # trsct_dir <- file.path(LOP_dir,'LOP_transects', fsep = .Platform$file.sep)
  # lidar.transect <- readOGR(dsn = trsct_dir, layer = "CAN_trsct_NAD1983LCC")
  # transect.buffer <- gBuffer(lidar.transect, width=400)
  # set.seed(paramsGL$global.seed)  ## spssample has a random component in where it places the grid, so needs to be initialized the same each time
  # val.hex <- spsample(transect.buffer, type = "hexagonal", cellsize = 50000)    ## cell size is set to 50km
  # 
  # val.hex.ID <- as.data.frame(seq(1, length(val.hex), by=1))
  # writeOGR( 
  #   
  #   SpatialPolygonsDataFrame( SpatialPolygons(val.hex), data=val.hex.ID, proj4string=CRS(proj4string(lidar.transect)) )
  #   
  #   
  #   , base_results_dir, 'val_hex', driver="ESRI Shapefile", overwrite_layer=TRUE)  ## just to check validity of other layers, never used later on
  # 
  # 
  # results.shp.df.reproj <- spTransform(results.shp.df, CRS(proj4string(val.hex)))
  # aaa <- over(results.shp.df.reproj, val.hex, returnList=FALSE, fn =mean)
  # 
  # ## generate hexagon unique IDs
  # Hex250ID <- as.data.frame(seq(1, length(HexPts.250m), by=1))
  # names(Hex250ID) <- "Hex250ID"  ## assign name to the column of the data frame
  # HexPts.250m.SPDF <- SpatialPointsDataFrame(coordinates(HexPts.250m), data=Hex250ID, proj4string=CRS(proj4string(lidar.2)))
  # 
  # #
  # #     # ggplot() + geom_point(data=temp.df, aes(x=Long, y=Lat, color=residual))
  # #
  # #     coordinates(temp.df) <- ~Long+Lat
  # #     proj4string(temp.df) <- CRS("+proj=longlat +datum=NAD83")
  ## TOO COMPLICATED
  
  
  ## round data and add nr of samples
  CAN.stats.df.rounded <- round(CAN.stats.df, 2)
  CAN.stats.df.rounded[, c("rsq", "ac_uns", "ac_sys")] <- round(CAN.stats.df[, c("rsq", "ac_uns", "ac_sys")], 3)
  CAN.stats.df.rounded[, "rmsdpct"] <- round(CAN.stats.df[, "rmsdpct"], 1)
  
  ECO.sorting.idx <- match(colnames(ECO.stats.df), as.character(names(TRN.samples.per.Ecozone))) ## to match order of list of Ecozones
  ECO.stats.df <- rbind(round(ECO.stats.df, 3), as.integer(TRN.samples.per.Ecozone)[ECO.sorting.idx], VAL.samples.per.Ecozone[ECO.sorting.idx])
  rownames(ECO.stats.df)[length(params3a$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
  
  UTM.sorting.idx <- match(as.character(zone.nrs$UTMzone), as.character(names(TRN.samples.per.UTMzone)))  ## to match order of list of UTM zones
  UTM.stats.df <- rbind(round(UTM.stats.df, 3), TRN.samples.per.UTMzone[UTM.sorting.idx], VAL.samples.per.UTMzone[UTM.sorting.idx])
  rownames(UTM.stats.df)[length(params3a$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
  
  CNG.sorting.idx <- match(params3a$Ch_attr.labels, as.character(names(TRN.samples.per.Ch_attr)))
  CNG.stats.df <- rbind(round(CNG.stats.df, 3), TRN.samples.per.Ch_attr[CNG.sorting.idx], VAL.samples.per.Ch_attr[CNG.sorting.idx])
  rownames(CNG.stats.df)[length(params3a$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
  
  ## stack results in df stats3a
  cmd <- sprintf('stats3a$%s$CAN.stats.df <- CAN.stats.df.rounded', method)
  eval(parse(text=cmd))
  cmd <- sprintf('stats3a$%s$ECO.stats.df <- ECO.stats.df', method)
  eval(parse(text=cmd))
  cmd <- sprintf('stats3a$%s$UTM.stats.df <- UTM.stats.df', method)
  eval(parse(text=cmd))
  cmd <- sprintf('stats3a$%s$CNG.stats.df <- CNG.stats.df', method)
  eval(parse(text=cmd))
  
  ## save csv and tex files
  cmd <- sprintf( "write.csv( stats3a$%s$CAN.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_CAN_stats.csv\", sep = '') )", method, method )
  eval(parse(text=cmd))
  cmd <- sprintf( "print( xtable(stats3a$%s$CAN.stats.df, digits=c(0,2,2,2,2,2,2,2,2,3,2,1,3,3,2)), file = file.path(base_results_dir, \"RES_MA_%s_CAN_stats.tex\", sep = '') )", method, method )
  eval(parse(text=cmd))
  
  cmd <- sprintf( "write.csv( stats3a$%s$ECO.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_ECO_%s.csv\", sep = '') )", method, method, params3a$metrics.MS )
  eval(parse(text=cmd))
  cmd <- sprintf( "print( xtable(stats3a$%s$ECO.stats.df, digits=c(0,rep(3, length(params3a$sampled.ecozones)))), file = file.path(base_results_dir, \"RES_MA_%s_ECO_%s.tex\", sep = '') )", method, method, params3a$metrics.MS)
  eval(parse(text=cmd))
  
  cmd <- sprintf( "write.csv( stats3a$%s$UTM.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_UTM_%s.csv\", sep = '') )", method, method, params3a$metrics.MS )
  eval(parse(text=cmd))
  cmd <- sprintf( "print( xtable(stats3a$%s$UTM.stats.df, digits=c(0,rep(3, length(paramsGL$zones)))), file = file.path(base_results_dir, \"RES_MA_%s_UTM_%s.tex\", sep = '') )", method, method, params3a$metrics.MS)
  eval(parse(text=cmd))
  
  cmd <- sprintf( "write.csv( stats3a$%s$CNG.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_CNG_%s.csv\", sep = '') )", method, method, params3a$metrics.MS )
  eval(parse(text=cmd))
  cmd <- sprintf( "print( xtable(stats3a$%s$CNG.stats.df, digits=c(0,rep(3, 4))), file = file.path(base_results_dir, \"RES_MA_%s_CNG_%s.tex\", sep = '') )", method, method, params3a$metrics.MS)
  eval(parse(text=cmd))
  
  ## Rsq barplots
  targs.for.barplots <- c("elev_p95", "cover_2m", "ag_biomass")
  
  ## by Ecozone
  data.to.plot <- t(ECO.stats.df[targs.for.barplots, ])
  data.melted <- melt(data.to.plot)
  colnames(data.melted) <- c("Ecozone", "Response_var", "R2")
  fig.name.str <- file.path(Assess.CAN.subdir, sprintf("Barplot_ECO_%s_%s.pdf", method, params3a$metrics.MS), sep='')
  pdf(fig.name.str)   ## , width=6, height=4
    theme_set(theme_gray(base_size = 18))
    barplot <- ggplot(data.melted, aes(Ecozone, R2)) +   
               geom_bar(aes(fill = Response_var), width = 0.75, position = "dodge", stat="identity") +
               theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                     axis.line = element_line(color="gray", size = 0.5),
                     panel.grid.major = element_line(colour="black", size=0.5, linetype="dashed"),
                     panel.grid.major.x = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_rect(colour = "gray", fill=NA),
                     legend.position="none") +
                     scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8)) +
               ylab("R^2")
  print(barplot)
  dev.off()
  
  
  ## by UTM zone
  data.to.plot <- t(UTM.stats.df[targs.for.barplots, ])
  data.melted <- melt(data.to.plot)
  colnames(data.melted) <- c("UTMzone", "Response_var", "R2")
  fig.name.str <- file.path(Assess.CAN.subdir, sprintf("Barplot_UTM_%s_%s.pdf", method, params3a$metrics.MS), sep='')
  pdf(fig.name.str, width=8, height=3)  
    theme_set(theme_gray(base_size = 18))
    barplot <- ggplot(data.melted, aes(UTMzone, R2)) +   
                geom_bar(aes(fill = Response_var), width = 0.75, position = "dodge", stat="identity") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                      axis.line = element_line(color="gray", size = 0.5),
                      panel.grid.major = element_line(colour="black", size=0.5, linetype="dashed"),
                      panel.grid.major.x = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "gray", fill=NA)) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8)) +
                ylab("R^2")
  print(barplot)
  dev.off()
  
  ## by Change attribution
  data.to.plot <- t(CNG.stats.df[targs.for.barplots, ])
  data.melted <- melt(data.to.plot)
  colnames(data.melted) <- c("ChangeAttribution", "Response_var", "R2")
  fig.name.str <- file.path(Assess.CAN.subdir, sprintf("Barplot_CNG_%s_%s.pdf", method, params3a$metrics.MS), sep='')
  pdf(fig.name.str)
    theme_set(theme_gray(base_size = 18))
    barplot <- ggplot(data.melted, aes(ChangeAttribution, R2)) +   
      geom_bar(aes(fill = Response_var), width = 0.75, position = "dodge", stat="identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.line = element_line(color="gray", size = 0.5),
            panel.grid.major = element_line(colour="black", size=0.5, linetype="dashed"),
            panel.grid.major.x = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "gray", fill=NA)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8)) +
      ylab("R^2")
  print(barplot)
  dev.off()

}

stats.file = file.path(base_results_dir, 'stats3a.Rdata', fsep = .Platform$file.sep) 
save(stats3a, file = stats.file)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
print(paste("Prog3a, total elapsed time:",seconds_to_period(toc[3])))