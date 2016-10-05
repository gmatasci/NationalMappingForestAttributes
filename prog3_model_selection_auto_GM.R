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
# - params3$run.MS <- T     
# - params3$run.MA <- T
# - params3$run.SG <- T 
# - params3$parallel.RF.MS <- T 
# - params3$parallel.RF.MA <- T (or F if we want RF importance)
# - params3$subsetting <- F
# - use with the right variables to exclude a priori

# - subset scatterplots to have lower density of points (colors convey that info) using something like: df <- df[df$dens<quantile(df$dens, 0.25),]


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
# -V yaImpute uses a RF with a nrTrees shared across all Ys, so actual nrTrees is nrTrees/nrYs -- solved by setting ntree=params3$ntree*length(params3$targ.yaImp.names.lg) in YaImp parameters
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
# -V reorder response variables to match order we want in paper
# -V run with final data: comparison between predicting response variables one at a time and not as a single multivariate Y: this means comparing yaimpute rf, rf, etc. -- done, we chose RF because of better performances ~ 0.1 R^2 higher than YaImpute
# -V remove from training set aberrant samples wrt each response variable (these samples will not impact the predictions, producing values in a meaningful range --> RF, Yaimp both respect the range of training set) or put in table 99.99% centile values -- not done, we'll keep it like this so far
# -V if parallel rf is used setting scale=t will have no effect, so use sequential rf for variable importance plot -- we use sequential RF, it takes a couple of hours to run with 100 trees
# -V procedure to see sampled values based on "FCID" or "POLY250ID" from X or Y dataframes loaded here linked to field "POLY250ID" of <UTMzone>_pt_centerpt_training_validation.shp (saved in /BAP_Imputation_working/wkg/<UTMzone>/)
#   then related to <UTMzone>_plot_inventory_attributes2.csv by "unique_id" -- this was the old way to do with Harold's FCIDs, now it is easier bc these IDs are the initial unique plot IDs of the LiDAR transect 
# -V if we only map on boreal do we exclude some ecozones? Shall we exclude them also in training and/or validation? -- if we use Brandt boreal definition, all of the 9 ecozones we sample are boreal ecozones.
# -V change bin width in residual hists -- set to 1m width for elev_p95 starting at 0

#### READS/WRITES ------------------------------------------------------------

## READS:
# - "lidar_metrics_mean_training_validation.csv": (from prog1) csv table with all observed LiDAR and forest attributes (Y) for selected TRN and VAL samples (3x3 polygons with average plot values or just single plot value) for all UTM zones
# - "poly_training_validation_exvars_extract.csv":  (from prog2) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

## WRITES:
# - "RF_<target>.Rdata": R data file with RF models for each response variable to use later for mapping (prediction) on the grid
# - "YAI.Rdata", "Ytrn.Rdata": R data file with YaImpute models and training set to use later for mapping (imputation) on the grid

#### INIT --------------------------------------------------------------------

print('Prog3: descriptive stats/plots, model selection and model assessment') 

rm(list=ls())

param_file <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatImp.R")

## subdirectory to save model assessment plots
Assess.CAN.subdir <- file.path(base_figures_dir, "Assessment_CAN_level", fsep = .Platform$file.sep)
if (! file.exists(Assess.CAN.subdir)){dir.create(Assess.CAN.subdir, showWarnings = F, recursive = T)}

#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params3 <- list()

## Actual parameters to be used
params3$subsetting <- T  ## to subset the dataset to a nr of samples = params3$nr.pts.plot (for debugging in development phase)
params3$methods <- list("YAI") ## accepts "RF" and/or "YAI", used to run analyses only for the specified methods
params3$run.checks <- F   ## whether to run checks for duplicate FCIDs, run only once, as of 06/06/2016 all is OK
params3$run.descr.stats <- F  ## whether to run descriptive stats block
params3$run.MS <- F     ## whether to run model selection block
params3$run.MA <- T    ## whether to run model assessment block (if set to FALSE, the script loads the prediction files saved in the last run in base_results_dir, to be used to change the plots based on the same results though)
params3$run.SG <- T    ## whether to run stats & graphs block
params3$parallel.RF.MS <- T    ## whether to run RF in parallel in the model selection block
params3$parallel.RF.MA <- F     ## whether to run RF in parallel in the model assessment block (variable importance is not available after having run the RF in parallel, so is set to FALSE) 
params3$ntree.MS <- 100     ## nr of RF trees in the model selection phase (ideally 500 but then YaImp, which cannot be run in parallel) it takes ages
params3$ntree.MA <- 100     ## nr of RF trees in the model assessment phase (very marginal improvements going higher than that, will cause bigger models and less priority on Westgrid for the mapping phase) 
params3$LongLat.val <- F   ## wheter to run the validation along longitudinal or latitudinal chunks of the transect 
params3$LongLat.val.nr.chunks <- 10   ##  number of chunks for long/lat validation

params3$mtry <- 'sqrt_nr_var'   ## how to set mtry parameter in RF
params3$nodesize <- 5       ## nodesize parameter of RF (5 is default for regression)
params3$plot.importance <- T   ## wheter to plot RF variable importance values

params3$predictor.groups <- list("Rthresh_0p95", "Rthresh_0p95_no_coords")   ## in model selection phase: list of predictors to include in the datasets for training and prediction
# params3$predictor.groups <- list("all", "SKstudy", "Rthresh_0p8", "Rthresh_0p95", "Rthresh_0p9_no_coords", "Rthresh_0p9_no_change", "Rthresh_0p9")  ## those used in a preceding exploratory analysis
params3$best.predictor.group <- "Rthresh_0p95"   ## in model assessment phase: final list of predictors
params3$predictors.to.rem.prior.know <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7', 'NDVI', 'Lat', 'Long')  ## list of predictors to remove a priori: all of the raw spectral bands (only removing VIS bands here results in the same set as b4, b5 and b7 are then filtered out based on |R|) and the NDVI (saturation issues)

params3$metrics.MS <- "rsq"    ## metrics that function regr_metrics() should return in model selection phase
params3$metrics.MA <- c("xmean", "xrangemin", "xrangemax", "xstdev", "ymean", "yrangemin", "yrangemax", "ystdev", "rsq", "rmsd", "rmsdpct", "ac_uns", "ac_sys", "bias")  ## metrics that function regr_metrics() should return in model assessment phase
params3$metrics.MA.colnames <- c("o_mean", "o_rangemin", "o_rangemax", "o_stdev", "p_mean", "p_rangemin", "p_rangemax", "p_stdev", "rsq", "rmsd", "rmsdpct", "ac_uns", "ac_sys", "bias")  ## relative column names for Latex tables

## final target variables to be predicted by both RF and YaImp
params3$targ.names.lg <- list("elev_mean", "elev_stddev", "elev_cv", "elev_p95", "percentage_first_returns_above_2m", "percentage_first_returns_above_mean", "loreys_height", "basal_area", "gross_stem_volume", "total_biomass")
params3$targ.names.sh <- list("el_m", "el_std", "el_cv", "el_p95", "p_1r_2m","p_1r_me", "loreys_h", "basal_a", "stem_v", "tot_biom")

## target variables given as input to YaImp to find best neighbor (other 4 are attached once FCIDs are known) 
params3$targ.yaImp.names.lg <- list("elev_mean", "elev_stddev", "elev_cv", "elev_p95", "percentage_first_returns_above_2m", "percentage_first_returns_above_mean")
params3$targ.yaImp.names.sh <- list("el_m", "el_std", "el_cv", "el_p95", "p_1r_2m", "p_1r_me")

## two key target variables used in descriptive plots to show X to Y relationship
params3$targ.key.names.lg <- list("percentage_first_returns_above_2m", "total_biomass")
params3$targ.key.names.sh <- list("p_1r_2m", "tot_biom")

## final target variable names and units to be used in plots
params3$targ.names.plots <- list("elev_mean", "elev_sd", "elev_cv", "elev_p95", "cover_2m", "cover_mean", "loreys_height", "basal_area", "stem_volume", "ag_biomass")
params3$targ.units.plots <- list("m", "m", "", "m", "%", "%", "m", "m^2/ha", "m^3/ha", "t/ha")

## retained 5 classes: 0="No change", 1="Fire", 2="Harvesting", 3="Non-stand replacing", 4="Infrastructure"
params3$Ch_attr.classes <- c("No change", "Fire", "Harvesting", "Non-stand replacing", "Infrastructure")
params3$Ch_attr.labels <- c(0, 1, 2, 3, 4)  ## no "Infrastructure" bc there are too few points

params3$groups.to.plot <- list('Bands', 'TCcomps_VI', 'ChangeAttribution', 'Years_Topo_Coords')  ## to plot relation with key response vars separately (and in a visually pleasing way)

params3$nr.pts.plot <- 2000    ## to make pairs plots less dense in points and to subset the dataset for tests in the developement phase

params3$corr.high <- 0.95    ## corr threshold to highlight high correlation among predictors
params3$corr.low <- 1     ## corr threshold to highlight low correlation between predictors and response vars

params3$nr.of.bins.UTMz <- 12  ## nr of bins for UTM zone histograms
params3$nr.of.bins.ecoz <- 18
params3$hist.centile <- 0.999  ## percentile threshold to subset data for histograms (avoid extreme values setting a too high xlim)

params3$distort.thresh.elevp95 <- 1.8     ## height threshold to show height over-under estimation wrt cover
params3$distort.thresh.cover2m <- 10    ## cover threshold to show height over-under estimation wrt height

## all ecozones in Canada
params3$all.ecozones <- list("Arctic Cordillera", "Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", 
                          "Hudson Plains", "Mixedwood Plains", "Montane Cordillera", "Northern Arctic", "Pacific Maritime", "Semiarid Prairies", 
                          "Southern Arctic", "Subhumid Prairies", "Taiga Cordillera", "Taiga Plains", "Taiga Shield East", "Taiga Shield West")

## ecozones covered by the transect
params3$sampled.ecozones <- list("Atlantic Maritime", "Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", "Hudson Plains", 
                                  "Taiga Plains", "Taiga Shield East", "Taiga Shield West")  ## Semiarid Prairies and Subhumid Prairies (ECOZONE_ID = 10 in shp) result in 0 and 3 plots, respectively, so we remove them

param_file_prog3 = file.path(base_wkg_dir, 'AllUTMzones_params3.Rdata', fsep = .Platform$file.sep) 
save(params3, file = param_file_prog3)

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

#### READ DATA --------------------------------------------------------------

## load plot training and validation data for lidar metrics and derived structural attributes (old "str", for structural)
Y.trn.val.raw <- read.csv(file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", fsep = .Platform$file.sep))

## change units for "summable" variables: unit per plot (625m^2) to unit per hectare (10000m^2) --> multiply by 16
Y.trn.val.raw[, c("basal_area", "gross_stem_volume", "total_biomass")] <- Y.trn.val.raw[, c("basal_area", "gross_stem_volume", "total_biomass")]*16

## transform biomass from kg/ha to t/ha
Y.trn.val.raw[, "total_biomass"] <- Y.trn.val.raw[, "total_biomass"]/1000

## load explanatory variables extracted for training plots (old "env", for environment)
X.trn.val.raw <- read.csv(file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", fsep = .Platform$file.sep))

X.trn.val.raw$Ch_attr <- as.factor(X.trn.val.raw$Ch_attr)   ## to make sure RF interprets Ch_attr as a categorical variable

## to check integrity of datasets before running everything
if (params3$run.checks) {
  
  ## check for order of sample points in the two datasets X and Y
  if ( !identical(X.trn.val.raw$FCID, Y.trn.val.raw$FCID) ) {
    stop("Prog3: order of FCIDs is different in lidar_metrics_mean_training_validation.csv wrt poly_training_validation_exvars_extract.csv")
  } 
  
  ## check for duplicate FCIDs in the dataset
  if ( length(unique(X.trn.val.raw$FCID)) != length(X.trn.val.raw$FCID) ) {
    stop("Prog3: FCIDs are not unique in lidar_metrics_mean_training_validation.csv and poly_training_validation_exvars_extract.csv")
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
    stop("Prog3: There are duplicate FCIDs (POLY250ID) within the various <UTMzone>_poly_training_validation.shp shapefiles")
  }
  if ( length(unique(FCIDs.df)) != length(FCIDs.df) ) {
    stop("Prog3: There are duplicate FCIDs (POLY250ID) between the various <UTMzone>_poly_training_validation.shp shapefiles")
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
stats3 <- list()
stats3$NA.samples.removed <- sum(idx.raster.zero==T) + sum(idx.Subhumid.Prairies==T)
TRN.samples.per.Ecozone <- table(Y.trn.val[Y.trn.val$TV=="TRAINING",]$ECOZONE)  ## infos to be added later on to the final stats dataframes
VAL.samples.per.Ecozone <- table(Y.trn.val[Y.trn.val$TV=="VALIDATION",]$ECOZONE)
TRN.samples.per.UTMzone <- table(X.trn.val[X.trn.val$TV=="TRAINING",]$UTMzone)
VAL.samples.per.UTMzone <- table(X.trn.val[X.trn.val$TV=="VALIDATION",]$UTMzone)
TRN.samples.per.Ch_attr <- table(X.trn.val[X.trn.val$TV=="TRAINING",]$Ch_attr)
VAL.samples.per.Ch_attr <- table(X.trn.val[X.trn.val$TV=="VALIDATION",]$Ch_attr)

#### CORRELATION MATRIX ----------------------------------------------

## compute the following outside of if block bc we want to know redundant predictors in all cases
prior.know.cont.idx <- !unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr) ## remove unwanted predictors by prior knowledge and not-continuous predictors
X.trn.val.corr <- X.trn.val[, unlist(params2$pred.names.sh)[prior.know.cont.idx]]
X.corr.matrix <- cor(X.trn.val.corr)  ## compute corr only among predictors (top-left sub-block)

## print correlation matrix as a Latex table
print( xtable(X.corr.matrix, digits=2), file = file.path(base_results_dir, "X_corr_matrix.tex", sep = '') )

## based on the dataset after removal of unwanted predictors by prior knowledge, find subset of redundant variables automatically:
## in a pair having r>cutoff remove the variable with highest average absolute correlation with rest of available variables (exact=T)
stats3$redund.pred.names.p80 <- findCorrelation(X.corr.matrix, cutoff=0.80, exact=T, verbose=F, names=T)
stats3$redund.pred.names.p90 <- findCorrelation(X.corr.matrix, cutoff=0.90, exact=T, verbose=F, names=T)
stats3$redund.pred.names.p95 <- findCorrelation(X.corr.matrix, cutoff=0.95, exact=T, verbose=F, names=T)

## descriptive statistics block
if (params3$run.descr.stats) {
  
  Y.trn.val.corr <- Y.trn.val[, unlist(params3$targ.names.lg)]   ## keep only response variables of interest
  colnames(Y.trn.val.corr) <- params3$targ.names.sh   ## rename them with short names
  
  X.corr.matrix.high <- X.corr.matrix
  X.corr.matrix.high[abs(X.corr.matrix) < params3$corr.high] <- NA  ## keep only high correlations among predictors
  
  data.matrix <- cbind(X.trn.val.corr, Y.trn.val.corr)   ## combined dataset including both Xs and Ys
  corr.matrix <- cor(data.matrix)
  XtoY.corr.matrix <- corr.matrix[1:ncol(X.trn.val.corr), seq(ncol(X.trn.val.corr)+1,ncol(corr.matrix),1)]  ## subset to keep only top-right corr matrix sub-block  
  XtoY.corr.matrix.low <- XtoY.corr.matrix
  XtoY.corr.matrix.low[abs(XtoY.corr.matrix) > params3$corr.low] <- NA   ## keep only low predictor-response correlations
  
  XY.corr.matrix <- round(cbind(X.corr.matrix, XtoY.corr.matrix), digits=3)  ## put together the two matrix sub-blocks
  XY.corr.matrix.NA <- round(cbind(X.corr.matrix.high, XtoY.corr.matrix.low), digits=3)
  
  ## save both complete matrix and those with NAs highlighting high correlations only
  write.csv(XY.corr.matrix, file = file.path(base_results_dir, "XY_corr_matrix.csv", sep = ''))
  str <- sprintf( "XY_corr_matrix_filtr_High%s_Low%s.csv", gsub(".", "p", as.character(params3$corr.high), fixed=T), gsub(".", "p", as.character(params3$corr.low), fixed=T)) 
  write.csv(XY.corr.matrix.NA, file = file.path(base_results_dir, str, sep = ''))

#### CANADIAN LEVEL X-Y RELATION PLOTS ------------------------------------------
  
  ## create figure subdirectories
  CAN.subdir <- file.path(base_figures_dir, "Relations_CAN_level", fsep = .Platform$file.sep)
  EcoUTMzone.subdir <- file.path(base_figures_dir, "Histograms_EcoUTMzone_level", fsep = .Platform$file.sep)
  if (! file.exists(CAN.subdir)){dir.create(CAN.subdir, showWarnings = F, recursive = T)}
  if (! file.exists(EcoUTMzone.subdir)){dir.create(EcoUTMzone.subdir, showWarnings = F, recursive = T)}
  
  ## create indices for subsampling 
  set.seed(paramsGL$global.seed)
  plot.idx.init <- sample(1:nrow(X.trn.val), params3$nr.pts.plot)  # take a subset of # params3$nr.pts.plot of points to plot
  UTM.idx <- which(paste('UTM', X.trn.val$UTMzone, sep='') %in% paramsGL$zones)  ## to only select samples of the UTM zones specified in paramsGL$zones
  plot.idx <- plot.idx.init[plot.idx.init %in% UTM.idx]
  
  Y.data.to.plot <- Y.trn.val[plot.idx, unlist(params3$targ.key.names.lg)] ## select Y columns to plot
  colnames(Y.data.to.plot) <- params3$targ.key.names.sh
  
  gr.idx <- 1
  for (pred.gr.names in params3$groups.to.plot) {
    gr.name <- unlist(params3$groups.to.plot[gr.idx])
    gr.idx <- gr.idx+1
    
    if (gr.name == 'ChangeAttribution') {   ## special case for ChangeAttribution
      pred.names <- "Ch_attr"
      X.data.to.plot <- as.data.frame(X.trn.val[, pred.names])  ## subset to keep only variable of interest
      colnames(X.data.to.plot) <- pred.names
      
      Y.data.to.plot.Ch_attr <- Y.trn.val[, unlist(params3$targ.key.names.lg)]
      colnames(Y.data.to.plot.Ch_attr) <- params3$targ.key.names.sh
      
      df.to.plot <- cbind(X.data.to.plot, Y.data.to.plot.Ch_attr)

      ## Canadian level (merged UTM zones) boxplots by Change_attr (categorical var)
      plot1 <- ggplot(df.to.plot, aes(x=Ch_attr, y=p_1r_2m, fill=Ch_attr)) +         ## canopy cover above 2m
        geom_boxplot(notch=F, outlier.shape = NA) + 
        scale_x_discrete(labels=params3$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plot2 <- ggplot(df.to.plot, aes(x=Ch_attr, y=tot_biom, fill=Ch_attr)) +          ## biomass
        geom_boxplot(notch=F, outlier.shape = NA) + 
        scale_x_discrete(labels=params3$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        coord_cartesian(ylim = c(0, 250000))
      str <- file.path(CAN.subdir, sprintf("CAN_boxplots_keyTarg_VS_%s.pdf", gr.name), sep='')
      pdf(str, width=10, height=5)
        grid.arrange(plot1, plot2, nrow=1, ncol=2, top=sprintf("%s VS pct_1ret_ab_2m and tot_biom", gr.name))  ## arrange the two plots in same two subplots
      dev.off()
      
      next
      
    } else if (gr.name == 'Bands') {  ## parameters for Bands
      pred.names <- unlist(params2$pred.names.sh$bands)
      cex.val <- 1    ## 
      cex.labels.val <- 1    ## size of variable labels
      cex.cor.val <- 0.8   ## size of text for correlation values
    } else if (gr.name == 'TCcomps_VI') {  ## parameters for Vegetation indices
      pred.names <- c(unlist(params2$pred.names.sh$TC), unlist(params2$pred.names.sh$VI))
      cex.val <- 1
      cex.labels.val <- 1
      cex.cor.val <- 0.8
    } else if (gr.name == 'Years_Topo_Coords') {   ## parameters for rest of variables
      pred.names <- c(unlist(params2$pred.names.sh$yrs), unlist(params2$pred.names.sh$topo), unlist(params2$pred.names.sh$coords))
      cex.val <- 1
      cex.labels.val <- 1
      cex.cor.val <- 0.8
    }
    
    ## correlation panels plots (common to all groups of features)
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
  

#### HISTOGRAMS BY UTM ZONE & ECOZONE ------------------------------------------
  
  ## continuous var
  cat.vars <- c("YrsSince_GrCh", "Ch_attr")  ## categorical vars to remove
  data.to.plot <- cbind(X.trn.val, Y.trn.val[, c(unlist(params3$targ.names.lg), "ECOZONE")])  ## not subset bc histograms can handle many points 

  for (var in c(unlist(params2$pred.names.sh), params3$targ.names.lg)) {   ## from column 6 onwards (5 is UTMzone and is the last one before the predictors) and till before-last bc last is ECOZONE
    
    if (!var %in% cat.vars) {  ## deal with continuous variables first...
      
      var.to.plot.compl <- data.to.plot[, var]
      ## indices of samples within some specified quantiles of the distribution (get rid of extremes)
      quant.idx <- var.to.plot.compl < quantile(var.to.plot.compl, params3$hist.centile, names=F) & var.to.plot.compl > quantile(var.to.plot.compl, 1-params3$hist.centile, names=F) ## index to select data within the 1-99th centile range to avoid extreme values skewing the histogram
      var.to.plot <- var.to.plot.compl[quant.idx] 
      UTMzone.to.plot <- data.to.plot[quant.idx, "UTMzone"]  ## subset of UTMzone and ECOZONE columns correspondig to those points
      ECOzone.to.plot <- data.to.plot[quant.idx, "ECOZONE"]
      my.hist.lims <- c(min(var.to.plot), max(var.to.plot))   ## set histogram limits and breaks that are the same across all the UTM and ecozones
      my.breaks <- c( seq(min(var.to.plot), max(var.to.plot), l=params3$nr.of.bins.UTMz+1) )
      
      ## by UTM zone
      str <- file.path(EcoUTMzone.subdir, sprintf("UTMzonesHist_%s.pdf", var), sep='')
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
      
      ## by Ecozone
      str <- file.path(EcoUTMzone.subdir, sprintf("EcozonesHist_%s.pdf", var), sep='')
      pdf(str)
      par(mfrow=c(3,3), oma = c(5,4,2,2) + 0.1, mar = c(2,2,2,2) + 0.1)  ## parameters to arrange nicely the panel of plots (for regular R plots, ggplot is different and requires package "gridExtra")
      for (z in 1:length(params3$sampled.ecozones)) {
        zone <- params3$sampled.ecozones[z]
        zone.idx <- ECOzone.to.plot == zone  ## get zone indexes to plot only data of that specific zone in this round of the loop
        hist(var.to.plot[zone.idx], main=sprintf("%s", zone), freq=T, breaks=my.breaks, xlim=my.hist.lims, xlab=NULL, ylab=NULL, col='lightblue')  # , breaks=20
      }
      title(var, outer=T)
      dev.off()
    
    } else {    ## ...then with categorical variables
      
      ## by UTM zone
      str <- file.path(EcoUTMzone.subdir, sprintf("UTMzonesHist_%s.pdf", var), sep='')
      pdf(str)
      par(mfrow=c(4,5), oma = c(5,4,2,2) + 0.1, mar = c(1,0,4.5,2) + 0.1)   
      for (z in 1:length(paramsGL$zones)) {
        zone <- paramsGL$zones[z]
        zone.nr <- substr(zone, 4, nchar(zone))
        zone.idx <- which(data.to.plot$UTMzone == zone.nr)
        if (var == "YrsSince_GrCh") {  ## special case for YrsSince_GrCh, some extra processing   
          yr.range <- 0:25  ## full range to be barplotted 
          yr.padded <- c( data.to.plot[zone.idx, var], yr.range[!yr.range %in% unique(data.to.plot[zone.idx, var])] )  ## add to the year data the years that are missing to have a complete set 
          barplot(table(yr.padded[yr.padded!=50]), main=sprintf("%s", zone), ylim=c(0, 300))  ## remove no-change data, years never appearing will have a value of 1 instead of 0 in the bar plot 
        } else if (var == "Ch_attr") {
          barplot(table(data.to.plot[zone.idx, var]), main=sprintf("%s", zone), ylim=c(0, 1200))
        }
      }
      title(var, outer=T)
      dev.off()
      
      ## by Ecozone
      str <- file.path(EcoUTMzone.subdir, sprintf("EcozonesHist_%s.pdf", var), sep='')
      pdf(str)
      par(mfrow=c(3,3), oma = c(5,4,2,2) + 0.1, mar = c(1,0,4.5,2) + 0.1)
      for (z in 1:length(params3$sampled.ecozones)) {
        zone <- params3$sampled.ecozones[z]
        zone.idx <- ECOzone.to.plot == zone
        if (var == "YrsSince_GrCh") {
          yr.range <- 0:25  ## full range to be barplotted
          yr.padded <- c( data.to.plot[zone.idx, var], yr.range[!yr.range %in% unique(data.to.plot[zone.idx, var])] )  ## add to the year data the years that are missing to have a complete set 
          barplot(table(yr.padded[yr.padded!=50]), main=sprintf("%s", zone), ylim=c(0, 300))  ## remove no-change data, years never appearing will have a value of 1 instead of 0 in the bar plot 
        } else if (var == "Ch_attr") {
          barplot(table(data.to.plot[zone.idx, var]), main=sprintf("%s", zone), ylim=c(0, 1200))
        }
      }
      title(var, outer=T)
      dev.off()
      
    }
    
  }

}

#### DATA PREPARATION -------------------------------------------------

## remove unneeded vars and predictors that we discard by a priori knowledge in X.trn.val data
stats3$initial.predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% params3$predictors.to.rem.prior.know]  ## save list of predictors at the start of the analysis (after removal of a priori exclusions) 
X.trn.val.data <- X.trn.val[, c('FCID', stats3$initial.predictors, 'TV')]

## extract core str response variables using in random forest analysis (leave out ECOZONE, to be retrieved later on through the FCID)
## imputation map only contains single band of plot FCIDs bc extracted key variables that get attached to map after imputation
Y.trn.val.data <- Y.trn.val[, c("FCID", unlist(params3$targ.names.lg), "TV")]

X.trn.val.data <- FCID_2_rownames(X.trn.val.data)  ## set FCID column as rownames and remove it from columns of df
Y.trn.val.data <- FCID_2_rownames(Y.trn.val.data)

## trick to be used in the testing phase only to have smaller datasets and a faster debugging
if (params3$subsetting) {
  set.seed(paramsGL$global.seed)
  plot.idx <- sample(1:nrow(X.trn.val.data), params3$nr.pts.plot)
  X.trn.val.data <- X.trn.val.data[plot.idx, ]
  Y.trn.val.data <- Y.trn.val.data[plot.idx, ]
}

## build the separate training and validation sets based on the TV column created in prog1
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

if (params3$run.MS) {  ## run model selection phase only if we want to select best combination of predictors that has then to be manually assigned to params3$best.predictor.group
  
  print('Model selection')
  
  ## loop to choose best subset of predictors
  RES_MS.df <- data.frame(matrix(nrow=0, ncol=length(params3$targ.names.lg)+1))  ## initialize empty df with as many columns as there are response variables (targets) plus one for the name of the group of predictors
  case.text <- sprintf("Method (ntree=%s, mtry=%s, nodesize=%s)", params3$ntree.MS, params3$mtry, params3$nodesize)  ## colname of first column is the details about which parameters of the RF are used (common to all)
  colnames(RES_MS.df) <- c(case.text, unlist(params3$targ.names.lg))

  for (predictor.gr in params3$predictor.groups) {
  
    print(sprintf('Predictors group: %s', predictor.gr))
    temp.tic <- proc.time() # start clocking time for each subset

    ## switch-case specifying different predictor combinaitons
    if (predictor.gr == "all") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% params3$predictors.to.rem.prior.know]  ## always remove the a priori excluded variables...
    } else if (predictor.gr == "SKstudy") {
      predictors <- c("TCB", "TCG", "TCW", "TCA", "TCD", "YrsSince_GrCh", "Ch_attr", "Elev", "Slope", "TWI", "TSRI", "Long", "Lat")  ## 13 predictors
    } else if (predictor.gr == "Rthresh_0p8") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p80)] ## ...plus those specified in stats3$redund.pred.names.pXX
    } else if (predictor.gr == "Rthresh_0p95_no_coords") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p95, params2$pred.names.sh$coords)]
    } else if (predictor.gr == "Rthresh_0p95") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p95)]
    } else if (predictor.gr == "Rthresh_0p9_no_coords") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p90, params2$pred.names.sh$coords)]
    } else if (predictor.gr == "Rthresh_0p9_no_change") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p90, params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr)]
    } else if (predictor.gr == "Rthresh_0p9") {
      predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p90)]
    }
    
    nr.vars <- length(predictors)
    
    ## set mtry in either of these two manners
    if (params3$mtry == 'sqrt_nr_var') { 
      mtries <- floor(sqrt(nr.vars))
    } else if (params3$mtry == 'nr_var_div_3') {
      mtries <- floor(nr.vars/3)
    }
  
    ## Random Forest (will be run anyways)
    RF.perf <- data.frame(matrix(nrow=1, ncol=length(params3$targ.names.lg)))  ## initialize row that will contain RF performances for each resp. variable
    colnames(RF.perf) <- params3$targ.names.lg
    
    if (params3$LongLat.val) { ## if we want to run the Long/Lat CV...
      targets.val <- "elev_p95"      ## ...we set the list of targets to iterate over with the next loop to only one variable (elev_p95)
    } else {
      targets.val <- params3$targ.names.lg  ## ..otherwise, in the standard case, it is the complete list
    }
    
    for (targ in targets.val) {
      
      if (!is.factor(X.trn$Ch_attr) | !is.factor(X.val$Ch_attr)) {  ## to double check it is interpreted as a categorical variable
        stop("Ch_attr is not set as a factor")
      }
      
      ## Training
      if (params3$parallel.RF.MS) {  ## in parallel with foreach()
        set.seed(paramsGL$global.seed)   ## set the same seed every time so that we obtain the same results
        nr.clusters <- min(params3$ntree.MS, detectCores())  ## set as minimum between the nr of trees and the nr of cores
        tot.nrtrees <- params3$ntree.MS      ## total nr of trees to be shared across cores (clusters)
        nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)   ## nr of trees per cluster (computed over the total nr of cluster minus 1)
        remaind <- tot.nrtrees%%(nr.clusters-1)   ## remainder of trees for last cluster
        cl <- makeCluster(nr.clusters)    ## initialize cores
        registerDoParallel(cl)   ## regirster them
        ## actual loop that will aggregate the results (output of all trees in the RF object rf.RF, the output of foreach()), use .multicombine=T for increased speed
        rf.RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages='randomForest') %dopar% {
          randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=ntrees, mtry=mtries, nodesize=params3$nodesize)
        }
        stopCluster(cl)
      } else {   ## or classical (sequential)
        set.seed(paramsGL$global.seed)
        rf.RF <- randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=params3$ntree.MS, mtry=mtries, nodesize=params3$nodesize)
      }
      
      ## Validation
      Y.val.predicted <- predict(rf.RF, X.val[,predictors], type="response", predict.all=F, nodes=F)
      if (!identical(rownames(X.val), rownames(Y.val))) {
        stop("TRN and VAL rownames do not match!!!")
      }
      ## apply function regr_metrics() on observed values Y.val$<targ> and predicted values Y.val.predicted 
      cmd=sprintf('perf <- as.vector(regr_metrics(Y.val$%s, as.matrix(Y.val.predicted))[[params3$metrics.MS]])', targ)   ## double [[]] to get only the value and not the name of list element
      eval(parse(text=cmd))
      RF.perf[1, targ] <- perf   ## fill RF.perf by column name targ
      
      ## Validation by CV on longitude/latitude splits (only for elev_p95)
      if (params3$LongLat.val && targ=="elev_p95") {
        
        LongLat.stats.elev_p95 <- data.frame( matrix( nrow=2, ncol=length(c(params3$metrics.MA, "b_yvsx", "a_yvsx")) ) )   ## initialize df to contain the stats (params3$metrics.MA) of CV...
        colnames(LongLat.stats.elev_p95) <- c(params3$metrics.MA.colnames, "b_yvsx", "a_yvsx")
        rownames(LongLat.stats.elev_p95) <- c("Long", "Lat")   ## ...by Long (1st row) and Lat (2nd row)
        for ( direction in c("Long", "Lat") ) {  ## outer loop to run the CV in the 2 directions independently
          
          ## initialize prediction vector to be filled at each CV round
          Y.val.predicted.LongLat <- data.frame(matrix(nrow=nrow(X.trn.val.data), ncol=1))
          rownames(Y.val.predicted.LongLat) <- rownames(X.trn.val.data)
          colnames(Y.val.predicted.LongLat) <- "elev_p95"
          
          ## set split var as either Long or Lat
          if (direction=="Long") {
            split.var <- X.trn.val.data$Long
          } else if (direction=="Lat") {
            split.var <- X.trn.val.data$Lat
          }
          
          ## compute params3$LongLat.val.nr.chunks quantiles 
          quantil <- quantile(split.var, seq(0, 1, 1/params3$LongLat.val.nr.chunks))
          for (i in 1:params3$LongLat.val.nr.chunks) {   ## loop over these chunks of transect
            if (i != params3$LongLat.val.nr.chunks) {   ## normal case where split.var value as to be within the 2 limits but without containing the upper
              idx.LongLat.val <- split.var >= quantil[i] & split.var < quantil[i+1]  ## set the logical indices for the left-out validation set for this CV round...
            } else {   ## special case to include max value in last split
              idx.LongLat.val <- split.var >= quantil[i] & split.var <= quantil[i+1]
            }
            idx.LongLat.trn <- !idx.LongLat.val  ## ...and set the indices for the kept-in training set as their negation
            rf.RF.LongLat <- randomForest(x=X.trn.val.data[idx.LongLat.trn, predictors], y=Y.trn.val.data[idx.LongLat.trn, targ], ntree=params3$ntree.MA, mtry=mtries, nodesize=params3$nodesize, importance=params3$plot.importance)
            prediction.res <- predict(rf.RF.LongLat, X.trn.val.data[idx.LongLat.val, predictors], type="response", predict.all=T, nodes=F)
            Y.val.predicted.LongLat[idx.LongLat.val,1] <- prediction.res$aggregate  ## retrieve actual RF prediction and assign it to the corresponding validation elements of the storing vector 
          }
          if ( !identical(rownames(Y.trn.val.data), rownames(Y.val.predicted.LongLat)) ) {
            stop("In LongLat validation, rownames for Y.trn.val.data and Y.val.predicted.LongLat are different!")
          }
          
          val.idx <- Y.trn.val.data$TV == "VALIDATION"  ## to compute the stats on the validation set only for better comparison with the classic validation approach (random 75/25% split over the transect)
          temp.df <- data.frame(observed=Y.trn.val.data$elev_p95[val.idx], predicted=Y.val.predicted.LongLat$elev_p95[val.idx])  ## temporary store the observed and predicted values for the validation set in the same df
          stat.row <- unlist(regr_metrics(temp.df$observed, temp.df$predicted)[c(params3$metrics.MA, "b_yvsx", "a_yvsx")])  ## compute assessment metrics and stats
          
          ## Long/Lat scatterplot of observed VS predicted values
          idx.targ <- which(params3$targ.names.plots==targ)
          fig.name.str <- file.path(Assess.CAN.subdir, sprintf("ObsPred_scatter_RF_%s_%s_%s.pdf", targ, direction, predictor.gr), sep='')
          title.str <- sprintf("%s [%s], %s split, %s", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ], direction, predictor.gr)
          box.str <- sprintf("R^2=%.3f, RMSE%%=%.1f, bias=%.2f", stat.row["rsq"], stat.row["rmsdpct"], stat.row["bias"])
          scatt.range <- c(min(temp.df$observed, temp.df$predicted), quantile(c(temp.df$observed, temp.df$predicted), 0.9999, names=F))   ## filter out the 99.99 percentile (outliers)
          GMFR.slope <- stat.row["b_yvsx"]
          GMFR.interc <- stat.row["a_yvsx"]
          pdf(fig.name.str)
          theme_set(theme_gray(base_size = 18))
          scatt.plot <- plot_colorByDensity(temp.df$observed, temp.df$predicted, xlim=scatt.range, ylim=scatt.range, xlab="observed", ylab="predicted", main=paste(title.str, '\n', box.str)) +
            geom_abline(intercept = GMFR.interc, slope = GMFR.slope, linetype="dashed", size=1) +   ## Add GMFR line
            geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line
          print(scatt.plot)
          dev.off()
          
          LongLat.stats.elev_p95[direction, ] <- stat.row   ## fill global stats df with results of this direction
          cmd=sprintf('stats3$LongLat.stats.elev_p95.%s <- LongLat.stats.elev_p95', predictor.gr)   ## save it in final list of results (stats3), specifying the predictor group used
          eval(parse(text=cmd))
        }
        
      }
      
    }
    
    RF.entry.name <- sprintf('RF.%s.nrvars%s', predictor.gr, nr.vars)  ## for the normal validation, save infos and put them as first element of RF.entry
    RF.entry <- as.list(c(RF.entry.name, RF.perf))
    list.RES_MS <- list(RES_MS.df, RF.entry)  ## save each RF.entry in a global list...
    RES_MS.df <- rbindlist(list.RES_MS)  ## ...to be binded by row to have the final df with the results
    
    ## YaImpute (to save on computing time will be run only for the selected group of predictors and if "YAI" is specified in the list of methods)
    if (predictor.gr == params3$best.predictor.group & "YAI" %in% params3$methods) {
      
      print('YaImpute')

      set.seed(paramsGL$global.seed) 
      
      ## - nodesize cannot be set
      ## - bootstrap=F to be used otherwise Error in yai(...: object 'xcvRefs' not found
      ## - rfMode="regression" to be used in all calls to yai(), otherwise it runs classification RF on discretized Ys (option "buildClasses") 
      ## - params3$ntree.MS*length(params3$targ.yaImp.names.lg) bc each Y uses a number of trees equal to ntree/nr.Ys
      yai.rf <- yai( x=X.trn[,predictors], y=Y.trn[,unlist(params3$targ.yaImp.names.lg)], bootstrap=F, method="randomForest", rfMode="regression", k=1, ntree=params3$ntree.MS*length(params3$targ.yaImp.names.lg), mtry=mtries)
      
      newtargets.output <- newtargets(yai.rf, X.val[,predictors], k=1)  ## run 1-NN (parameter k=1) imputation (analogue to prediction) on validation set...
      ids.val.predicted <- newtargets.output$neiIdsTrgs   ## ...and get the predicted IDs (FCIDs of the training samples)
      ids.val.predicted <- rownames_2_FCID(ids.val.predicted)  ## add column with IDs of the validation samples
      colnames(ids.val.predicted)[2] <- "Predicted_FCID"
      
      Y.trn.to.assign <- rownames_2_FCID(Y.trn)  ## add column of sample IDs to the training set
      
      ## join by FCID the two dataframes to assign imputed values from the training set
      Y.val.predicted <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y = "FCID")
      rm(Y.trn.to.assign)
      
      ## sort both observed and predicted dataframes to have same row order
      Y.val.predicted <- arrange(Y.val.predicted, FCID)
      Y.val.observ <- rownames_2_FCID(Y.val)
      Y.val.observ <- arrange(Y.val.observ, FCID)
      
      yaImp.perf <- vector('double')  ## initialize empty vector to store metrics
      for (targ in params3$targ.names.lg) {   ## loop over the resp. variables to compute the assessment metrics for the model selection phase (params3$metrics.MS) 
        cmd=sprintf('perf <- as.vector(regr_metrics(Y.val.observ$%s, as.matrix(Y.val.predicted$%s))[[params3$metrics.MS]])', targ, targ)   ## double [[]] to get only the value and not the name of list element
        eval(parse(text=cmd))
        yaImp.perf <- c(yaImp.perf, perf)
      }
      rm(Y.val.observ)
      
      ## fill final df with the results (same one as for RF)
      yaImp.entry.name <- sprintf('yaImp.%s.nrvars%s', params3$best.predictor.group, nr.vars)
      yaImp.entry <- as.list(c(yaImp.entry.name, yaImp.perf))
      list.RES_MS <- list(RES_MS.df, yaImp.entry)
      RES_MS.df <- rbindlist(list.RES_MS) 
      
    }

    temp.toc <- proc.time()-temp.tic[3]
    print(paste(predictor.gr, "elapsed time:", seconds_to_period(temp.toc[3])))
      
  }
  
  RES_MS.df <- t(RES_MS.df)  ## transpose df for better visualization and save it as csv file
  write.csv(RES_MS.df, file = file.path(base_results_dir, sprintf("RES_MS_%s.csv", params3$metrics.MS), sep = ''))

  stats3$RES_MS.df <- RES_MS.df  ## add table to stats list to be saved
  
}

#### MODEL ASSESSMENT ---------------------------------------------------------------------------

## create model subdirectories to save models for later use (in mapping phase)
models.subdir <- file.path(base_wkg_dir, "Models", fsep = .Platform$file.sep)
if (! file.exists(models.subdir)){dir.create(models.subdir, showWarnings = F, recursive = T)}

## specify best set of predictors found after model selction phase (params3$best.predictor.group): allow selection only among plausible models (e.g. Rthresh_0p8 discards Latitude)
if (params3$best.predictor.group == "all") {
  predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% params3$predictors.to.rem.prior.know]
} else if (params3$best.predictor.group == "Rthresh_0p9") {
  predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p90)]
} else if (params3$best.predictor.group == "Rthresh_0p95") {
  predictors <- unlist(params2$pred.names.sh)[!unlist(params2$pred.names.sh) %in% c(params3$predictors.to.rem.prior.know, stats3$redund.pred.names.p95)]
}

nr.vars <- length(predictors)

stats3$final.predictors <- predictors   ## save final set of predictors

if (params3$run.MA) {  ## run this block only if you want to (re)run the actual analysis, otherwise, to just redo plots, set params3$run.MA <- F and params3$run.SG <- T ... (see else block)  
  
  print('Model assessment')
  
  if (params3$mtry == 'sqrt_nr_var') {
    mtries <- floor(sqrt(nr.vars))
  } else if (params3$mtry == 'nr_var_div_3') {
    mtries <- floor(nr.vars/3)
  }

  for (method in params3$methods) {

#### RANDOM FOREST ---------------------------------------------------------------------
    
    if (method == "RF") {  ## run only if RF is in params3$methods
      
    print('Random Forest')
    temp.tic <- proc.time()    ## start specific clock for RF
    
      Y.val.predicted.RF <- data.frame(matrix(nrow=nrow(X.val), ncol=length(params3$targ.names.lg))) ## initialize df to store predictions: nr of rows equal to nr of validation samples, nr of column equal to nr of resp. variables
      rownames(Y.val.predicted.RF) <- rownames(X.val)
      colnames(Y.val.predicted.RF) <- params3$targ.names.lg
      Y.val.predicted.stddev.RF <- Y.val.predicted.RF  ## initialize the same way df to store standard deviations over the RF trees (measure of uncertainty)
      for (targ in params3$targ.names.lg) {
        
        if (!is.factor(X.trn$Ch_attr) | !is.factor(X.val$Ch_attr)) {
          stop("Ch_attr is not set as a factor")
        }
        
        ## Training (same as for MS phase but with params3$<param>.MA parameters)
        if (params3$parallel.RF.MA) {
          set.seed(paramsGL$global.seed)
          nr.clusters <- min(params3$ntree.MA, detectCores())
          tot.nrtrees <- params3$ntree.MA
          nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)
          remaind <- tot.nrtrees%%(nr.clusters-1)
          cl <- makeCluster(nr.clusters)
          registerDoParallel(cl)
          rf.RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages='randomForest') %dopar% {
            randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=ntrees, mtry=mtries, nodesize=params3$nodesize, importance=params3$plot.importance)
          }
          stopCluster(cl)
        } else {
          set.seed(paramsGL$global.seed)
          rf.RF <- randomForest(x=X.trn[,predictors], y=Y.trn[,targ], ntree=params3$ntree.MA, mtry=mtries, nodesize=params3$nodesize, importance=params3$plot.importance)
        }
        
        ## save RF model for each target in a R object (to be used in IDL scripts for the actual mapping)
        cmd <- sprintf('RF.model.path <- file.path(models.subdir, "RF_%s.Rdata", fsep=.Platform$file.sep)', targ)
        eval(parse(text=cmd))
        save(rf.RF, file=RF.model.path)
        
        ## validation
        prediction.res <- predict(rf.RF, X.val[,predictors], type="response", predict.all=T, nodes=F)
        Y.val.predicted.RF[,targ] <- prediction.res$aggregate  ## save actual predictions in the corresponding column for each resp. variable
        Y.val.predicted.stddev.RF[,targ] <- apply(prediction.res$individual, 1, sd)  ## save std dev over the RF trees (std dev of the individual prediction of each tree, a row of prediction.res$individual)
        
      }
      
      ## arrange and save predictions and uncertainties to be loaded in case we want only to redo plots
      Y.val.predicted.RF <- rownames_2_FCID(Y.val.predicted.RF)  ## to make the column FCID available for the merge() function in the assessment part below
      Y.val.predicted.RF <- arrange(Y.val.predicted.RF, FCID)
      Y.val.predicted.stddev.RF <- rownames_2_FCID(Y.val.predicted.stddev.RF)  ## to make the column FCID available for the merge() function in the assessment part below
      Y.val.predicted.stddev.RF <- arrange(Y.val.predicted.stddev.RF, FCID)
      write.csv(Y.val.predicted.RF, file = file.path(base_results_dir, sprintf("Y_val_predicted_matrix_RF.csv"), sep = ''))
      write.csv(Y.val.predicted.stddev.RF, file = file.path(base_results_dir, sprintf("Y_val_predicted_stddev_matrix_RF.csv"), sep = ''))
      
      ## compute processing time for RF
      temp.toc <- proc.time()-temp.tic[3]
      print(paste("Random Forest elapsed time:", seconds_to_period(temp.toc[3])))
    
    }
  

#### YAIMPUTE -------------------------------------------------------------------
    
    if (method == "YAI") {  ## run only if "YAI" is specified in methods
    
      print('YaImpute')
      temp.tic <- proc.time()
      
      set.seed(paramsGL$global.seed)
      ## - nodesize cannot be set
      ## - bootstrap=F to be used otherwise Error in yai(...: object 'xcvRefs' not found
      ## - rfMode="regression" to be used in all calls to yai(), otherwise it runs classification RF on discretized Ys (option "buildClasses") 
      ## - params3$ntree.MS*length(params3$targ.yaImp.names.lg) bc each Y uses a number of trees equal to ntree/nr.Ys
      yai.rf <- yai( x=X.trn[,predictors], y=Y.trn[,unlist(params3$targ.yaImp.names.lg)], bootstrap=F, method="randomForest", rfMode="regression", k=1, ntree=params3$ntree.MA*length(params3$targ.yaImp.names.lg), mtry=mtries)
      
      ## variable importance plot: commented because it is not clear how it is computed for imputation
      # Yai.var.imp <- yaiVarImp(yai.rf, nTop =0, plot=F)
      # par(las=1)
      # par(mar=c(4,15,1,1))
      #   boxplot(Yai.var.imp, horizontal = TRUE, xlab = "Scaled Variable Importance")
      # dev.off()

      ## impute IDs (same as for model selection phase)
      newtargets.output <- newtargets(yai.rf, X.val[,predictors], k=1)
      
      # ## ------- TODEL ---
      # ## TEST IF ANN PARAMETER HAS IMPACT
      # temp.toc <- proc.time()-temp.tic[3]
      # print(paste("usual elapsed time:", seconds_to_period(temp.toc[3])))
      # print(paste("usual Yai obj size:", object.size(yai.rf)))
      # print(paste("usual newtargets.output obj size:", object.size(newtargets.output)))
      # 
      # temp.tic <- proc.time()
      # set.seed(paramsGL$global.seed)
      # yai.rf.ANN <- yai( x=X.trn[,predictors], y=Y.trn[,unlist(params3$targ.yaImp.names.lg)], bootstrap=F, method="randomForest", rfMode="regression", k=1, ntree=params3$ntree.MA*length(params3$targ.yaImp.names.lg), mtry=mtries, ann=T)
      # newtargets.output.ANN <- newtargets(yai.rf.ANN, X.val[,predictors], k=1, ann=T)
      # ids.val.predicted.ANN <- newtargets.output.ANN$neiIdsTrgs
      # 
      # temp.toc <- proc.time()-temp.tic[3]
      # print(paste("with ANN elapsed time:", seconds_to_period(temp.toc[3])))
      # print(paste("with ANN Yai obj size:", object.size(yai.rf.ANN)))
      # print(paste("with ANN newtargets.output obj size:", object.size(newtargets.output.ANN)))
      # ## ------- TODEL ---
      
      ids.val.predicted <- newtargets.output$neiIdsTrgs
      ids.val.predicted <- rownames_2_FCID(ids.val.predicted)
      colnames(ids.val.predicted)[2] <- "Predicted_FCID"
      
      Y.trn.to.assign <- rownames_2_FCID(Y.trn)
      
      ## join by FCID the two dataframes to get imputed values
      Y.val.predicted.Yai <- merge(ids.val.predicted, Y.trn.to.assign, by.x="Predicted_FCID", by.y="FCID")
      Y.val.predicted.Yai$Predicted_FCID <- NULL
      
      ## save model (how to assign IDs) and Ytrn (from where to copy the values to paste) to impute values when mapping in IDL
      Yai.model.path <- file.path(models.subdir, "YAI.Rdata", fsep=.Platform$file.sep)   
      save(yai.rf, file=Yai.model.path)
      Y.trn.path <- file.path(models.subdir, "Ytrn.Rdata", fsep=.Platform$file.sep)
      save(Y.trn.to.assign, file=Y.trn.path)
      rm(Y.trn.to.assign)
      
      ## arrange and save predictions to be loaded in case we want only to redo plots
      Y.val.predicted.Yai <- arrange(Y.val.predicted.Yai, FCID)
      write.csv(Y.val.predicted.Yai, file = file.path(base_results_dir, sprintf("Y_val_predicted_matrix_YAI.csv"), sep = ''))
      
      ## compute processing time for YAI
      temp.toc <- proc.time()-temp.tic[3]
      print(paste("YaImpute elapsed time:", seconds_to_period(temp.toc[3])))
    
    }
  
  }  ## end for on params3$methods

} else {  ## ...if we did not run the model assessment, read the csv with the predictions that we saved beforehand
  
  if (params3$run.SG) {   ## only run this block if we want to run the STATS AND GRAPHS part below
    for (method in params3$methods) {
      if (method == "RF") {  ## only read the csv files of the corresponding methods
        Y.val.predicted.RF <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_matrix_RF.csv"), fsep = .Platform$file.sep), row.names=1)  ## row.names=1 to have same format as when these dataframes are produced in the predictions phase
        Y.val.predicted.stddev.RF <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_stddev_matrix_RF.csv"), fsep = .Platform$file.sep), row.names=1)
      }
      if (method == "YAI") {
        Y.val.predicted.Yai <- read.csv(file.path(base_results_dir, sprintf("Y_val_predicted_matrix_YAI.csv"), fsep = .Platform$file.sep), row.names=1)
      }
    }
  }
  
}  ## end if-else on params3$run.MA

#### STATS AND GRAPHS ---------------------------------------------------------------

if (params3$run.SG) {  ## only run this block if we want to run this STATS AND GRAPHS part

  ## if it doesn't already exist, create folder to save bivariate check plots
  BivCheck.CAN.subdir <- file.path(Assess.CAN.subdir, "Bivariate_check", sep='')  
  if (! file.exists(BivCheck.CAN.subdir)){dir.create(BivCheck.CAN.subdir, showWarnings = F, recursive = T)}
  
  ## same for residuals plots
  Resid.subdir <- file.path(Assess.CAN.subdir, "Residuals", sep='')  
  if (! file.exists(Resid.subdir)){dir.create(Resid.subdir, showWarnings = F, recursive = T)}
  
  for (method in params3$methods) {  ## again go over the list of methods...
  
    if (method == "RF") {
      Y.val.predicted <- Y.val.predicted.RF   ## ...and assign the corresponding df of predicted values to the general df "Y.val.predicted"
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
    CAN.stats.df <- data.frame( matrix( nrow=length(params3$targ.names.plots), ncol=length(params3$metrics.MA.colnames) ) ) ## initialize matrix to store global assessment stats
    colnames(CAN.stats.df) <- params3$metrics.MA.colnames  ## columns are the metrics
    rownames(CAN.stats.df) <- params3$targ.names.plots   ## rows are resp. variable names
    
    ECO.stats.df <- data.frame( matrix( nrow=length(params3$targ.names.plots), ncol=length(params3$sampled.ecozones) ) )  ## initialize matrix to store by ecozone assessment stats
    colnames(ECO.stats.df) <- params3$sampled.ecozones  ## columns are the sampled ecozones
    rownames(ECO.stats.df) <- params3$targ.names.plots  
    
    UTM.stats.df <- data.frame( matrix( nrow=length(params3$targ.names.plots), ncol=length(paramsGL$zones) ) )  ## initialize matrix to store by UTM zone assessment stats
    colnames(UTM.stats.df) <- paramsGL$zones    ## columns are the UTM zones
    rownames(UTM.stats.df) <- params3$targ.names.plots 
    
    CNG.stats.df <- data.frame( matrix( nrow=length(params3$targ.names.plots), ncol=length(params3$Ch_attr.classes)-1 ) )  ## initialize matrix to store by Change attribution assessment stats (we do not consider the class "Infrastructure")
    colnames(CNG.stats.df) <- params3$Ch_attr.classes[-length(params3$Ch_attr.classes)]
    rownames(CNG.stats.df) <- params3$targ.names.plots
    
    ## define axis labels for plots
    plot.xlabel <- "observed"
    plot.ylabel <- "predicted"
    plot.xlabel.resid <- "predicted"
    plot.ylabel.resid <- "residuals"
    
    ## initialize empty dataframe to host the complete results table to be saved in a shp: nrow = nr of validation samples, ncol = 10 targ * 3 values (observed, predicted, residual) = 30
    results.shp.df <- data.frame( matrix( nrow=nrow(Y.val.observ), ncol=3*length(params3$targ.names.lg) ) ) 
    rownames(results.shp.df) <- Y.val.observ$FCID
    
    ## initialize variable importance dataframe
    var.imp.df <- data.frame( matrix( nrow=length(params3$targ.names.lg), ncol=nr.vars ) ) ## initialize matrix to store variable importance values
    colnames(var.imp.df) <- predictors
    rownames(var.imp.df) <- params3$targ.names.lg
    

#### ACCURACY ASSESSMENT ---------------------------------------------------------------
    
    idx.targ <- 1  ## growing numerical index to automatically fill results.shp.df and other dfs
    for (targ in params3$targ.names.lg) {
      
      ## load RF model for each target to get variable importance
      cmd <- sprintf('load(file.path(models.subdir, "RF_%s.Rdata", fsep=.Platform$file.sep))', targ) 
      eval(parse(text=cmd))
  
      ## fill variable importance dataframe
      ## type=1: mean decrease in accuracy (mean increase in MSE for regression). For each tree, the prediction error on the out-of-bag portion of the data is
      ## recorded. Then the same is done after permuting each predictor variable. The difference between the two are then averaged over all trees, and normalized by the standard deviation of the differences.
      ## IF PARALLEL RF IS USED SETTING SCALE=T WILL HAVE NO EFFECT, SO USE SEQUENTIAL RF FOR VARIABLE IMPORTANCE PLOT
      var.imp.df[targ, ] <- t(importance(rf.RF, type=1, scale=T))  
      
      ## build temporary df to enable stats by category with ddply()
      cmd <- sprintf('temp.df <- data.frame(Y.val.observ$%s, Y.val.predicted$%s, Y.val.predicted[, c("ECOZONE", "UTMzone")], X.val[, c("Ch_attr", "Long", "Lat")])', targ, targ)
      eval(parse(text=cmd))
      colnames(temp.df) <- c("observed", "predicted", "ECOZONE", "UTMzone", "Ch_attr", "Long", "Lat")
      temp.df$resid <- temp.df$predicted - temp.df$observed  ## compute residuals
      
      ## store in results.shp.df all predictions and residuals to be written in a shp
      results.shp.df[, seq((idx.targ-1)*3+1, (idx.targ-1)*3+3)] <- data.frame(temp.df$observed, temp.df$predicted, temp.df$resid)  ## fill df 3 columns at a time
      cmd <- sprintf('colnames(results.shp.df)[seq((idx.targ-1)*3+1, (idx.targ-1)*3+3)] <- c(\"O_%s\", \"P_%s\", \"R_%s\")', params3$targ.names.sh[idx.targ],  params3$targ.names.sh[idx.targ],  params3$targ.names.sh[idx.targ])
      eval(parse(text=cmd))
  
      ## CAN-level stats
      CAN.stats.compl <- unlist(regr_metrics(temp.df$observed, temp.df$predicted)[c(params3$metrics.MA, "b_yvsx", "a_yvsx")])  ## specified metrics to report in summary table in the paper plus a and b coefficients for GMFR
      CAN.stats.df[idx.targ, ] <- as.vector(CAN.stats.compl[params3$metrics.MA])
        
      ## CAN-level scatterplots
      fig.name.str <- file.path(Assess.CAN.subdir, sprintf("ObsPred_scatter_%s_%s.pdf", method, targ), sep='')
      title.str <- sprintf("%s [%s]: ", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ])
      box.str <- sprintf("R^2=%.3f, RMSE%%=%.1f, bias=%.2f", CAN.stats.df[idx.targ, "rsq"], CAN.stats.df[idx.targ, "rmsdpct"], CAN.stats.df[idx.targ, "bias"])
      scatt.range <- c(min(temp.df$observed, temp.df$predicted), quantile(c(temp.df$observed, temp.df$predicted), 0.9999, names=F))
      GMFR.slope <- CAN.stats.compl["b_yvsx"]
      GMFR.interc <- CAN.stats.compl["a_yvsx"]
      ## uncomment to check if corr. coeff. R is equal to the slope of the regression based on standardized predicted and observed values: pred = a + b obs
      # lm.mod <- lm(predicted~observed, data=temp.df)
      # lm.slope <- lm.mod$coefficients[2]
      # lm.interc <- lm.mod$coefficients[1]
      # temp.df.stand <- as.data.frame(scale(temp.df[, c("predicted", "observed")], center=T, scale=T))
      # lm.mod.stand <- lm(predicted~observed, data=temp.df.stand)
      # if ( as.vector(lm.mod.stand$coefficients[2]) - as.vector(sqrt(CAN.stats.df[idx.targ, "rsq"])) > 1e-10 ) {
      #   warning(sprintf("Predicted VS observed R and standardized regression slope do not match for %s on %s", method, targ))
      # }
      pdf(fig.name.str)
        theme_set(theme_gray(base_size = 18))
        scatt.plot <- plot_colorByDensity(temp.df$observed, temp.df$predicted, xlim=scatt.range, ylim=scatt.range, xlab=plot.xlabel, ylab=plot.ylabel, main=paste(title.str, '\n', box.str)) +
                    geom_abline(intercept = GMFR.interc, slope = GMFR.slope, linetype="dashed", size=1) +   ## Add GMFR line
                    # geom_abline(intercept = lm.interc, slope = lm.slope, linetype="dashed", size=1, color='red') +   ## Add lm line
                    geom_abline(intercept = 0, slope = 1, colour = 'black', size=1)  # Add 1:1 line line
      print(scatt.plot)
      dev.off()
  
      ## by Ecozone stats
      stats.by <- ddply(temp.df, "ECOZONE", summarize, regr_metrics(observed, predicted)[params3$metrics.MS])   ## apply regr_metrics to the columns "predicted" and "observed" of dataframe temp.df split by "ECOZONE"...
      ECO.stats.df[idx.targ, ] <- t(stats.by[match(as.character(params3$sampled.ecozones), as.character(stats.by$ECOZONE)), ][,2])   ## ...then match the order of the stats.by dataframe to the list of Ecozones (in params3) and finally take only 2nd column (with results) and transpose it to fill ECO.stats.df 
      
      ## by UTMZone stats
      zone.nrs <- data.frame("UTMzone"=substr(paramsGL$zones, 4, nchar(paramsGL$zones)))  ## to get rid of the "UTM" characters
      stats.by <- ddply(temp.df, "UTMzone", summarize, regr_metrics(observed, predicted)[params3$metrics.MS])
      UTM.stats.df[idx.targ, ] <- t(stats.by[match(as.character(zone.nrs$UTMzone), as.character(stats.by$UTMzone)), ][,2])
      
      ## by Ch_attr stats
      temp.df.Ch <- temp.df[temp.df$Ch_attr!=4, ]   ## to drop "infrastructure" (level nr 4) 
      temp.df.Ch$Ch_attr <- factor(temp.df.Ch$Ch_attr)  ## to recompute factor level after having dropped one level
      stats.by <- ddply(temp.df.Ch, "Ch_attr", summarize, regr_metrics(observed, predicted)[params3$metrics.MS])   
      CNG.stats.df[idx.targ, ] <- t(stats.by[match(params3$Ch_attr.labels, stats.by$Ch_attr), ][,2])
      
      ## histograms of residuals by Ecozone
      if (targ %in% c("elev_p95", "elev_mean", "loreys_height")) {   ## if resp. variable is measured in meters start around 0 and set 1m wide bins
        idx.to.plot <- between(temp.df$resid, lower=quantile(temp.df$resid, 0.001), upper=ceiling(quantile(temp.df$resid, 0.999) ) )  ## remove outliers
        my.hist.lims <- c(floor(quantile(temp.df$resid, 0.001)), ceiling(quantile(temp.df$resid, 0.999)))  ## round to the lower integer the quantiles (min and max values)
        my.breaks <- seq(my.hist.lims[1], my.hist.lims[2], l=sum(abs(my.hist.lims))+1 )  ## set bin breaks at each integer between min and max
      } else {   
        idx.to.plot <- rep(T, length(temp.df$resid))
        my.hist.lims <- c(min(temp.df$resid), max(temp.df$resid))  
        my.breaks <- c( seq(min(temp.df$resid), max(temp.df$resid), l=params3$nr.of.bins.ecoz+1) )
      }
      str <- file.path(Resid.subdir, sprintf("ResidualsHist_%s_%s.pdf", method, targ), sep='')
      pdf(str)
      par(mfrow=c(3,3), oma = c(5,4,2,2) + 0.1, mar = c(2,2,2,2) + 0.1)  ## parameters to arrange nicely the panel of plots
      for (z in 1:length(params3$sampled.ecozones)) {
        zone <- params3$sampled.ecozones[z]
        zone.idx <- temp.df$ECOZONE == zone  ## get zone indexes to plot only data of that specific zone in this round of the loop
        hist(temp.df$resid[zone.idx & idx.to.plot], main=sprintf("%s", zone), freq=T, breaks=my.breaks, xlim=my.hist.lims, xlab=NULL, ylab=NULL, col='red')  # , breaks=20
      }
      title(targ, outer=T)
      dev.off()
      
      ## scatterplots of predicted values (now this is the x-axis) VS residuals at CAN level
      fig.name.str <- file.path(Resid.subdir, sprintf("ResidPred_scatter_%s_%s.pdf", method, targ), sep='')
      title.str <- sprintf("%s [%s]: ", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ])
      pdf(fig.name.str)
        theme_set(theme_gray(base_size = 18))
        resid.plot <- plot_colorByDensity(temp.df$predicted, temp.df$resid, xlim=scatt.range, xlab=plot.xlabel.resid, ylab=plot.ylabel.resid, main=title.str) +
          geom_abline(intercept = 0, slope = 0, colour = 'black', size=1, linetype="dashed")  # Add horizontal line line
      print(resid.plot)
      dev.off()
      
      idx.targ <- idx.targ+1
      
    }   ## end for targ in params3$targ.names.lg

#### PLOT VARIABLE IMPORTANCE ---------------------------------------------------------------
    
    if (params3$plot.importance) {  ## run only if specified in params3$plot.importance and only for RF
      
      rownames(var.imp.df) <- params3$targ.names.plots
      
      ## overall boxplot over the 10 targets
      fig.name.str <- file.path(Assess.CAN.subdir, sprintf("RF_VariableImportance_IncMSE_overall.pdf"), sep='')
      pdf(fig.name.str)
        mar.default <- c(5,4,4,2) + 0.1
        par(mar = mar.default + c(0, 4, 0, 0)) 
        theme_set(theme_gray(base_size = 18))
        medians <- sapply(var.imp.df, median)  ## compute median values
        idx.medians <- order(medians, decreasing=T)   ## sort by decreasing median value
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
    
    
#### DISTORTION OF COVARIANCE ---------------------------------------------------------------
    
    ## flag of critical samples with over-under estimation in two key variables (elev_p95 vs cover_2m)
    results.shp.df$OvHtUnCov <- results.shp.df$R_el_p95 > params3$distort.thresh.elevp95 & results.shp.df$R_p_1r_2m < -params3$distort.thresh.cover2m ## overestimation of height and underestimation of cover
    results.shp.df$UnHtOvCov <- results.shp.df$R_el_p95 < -params3$distort.thresh.elevp95 & results.shp.df$R_p_1r_2m > params3$distort.thresh.cover2m  ## underestimation of height and overestimation of cover
    
    ## write shp with individual validation plot prediction results
    results.shp.df <- data.frame(results.shp.df, X.val$Long, X.val$Lat, rownames(results.shp.df))   ## add FCID as a column of the dataframe to be able to read in back in ArcMap or in R with readOGR
    colnames(results.shp.df)[ncol(results.shp.df)] <- "FCID"   ## rename it to respect the 10 char limits in ESRI shapefiles!
    coordinates(results.shp.df) <- ~X.val.Long+X.val.Lat   ## set coordinates as center of the plot
    proj4string(results.shp.df) <- CRS("+proj=longlat +datum=NAD83")
    writeOGR(results.shp.df, base_results_dir, sprintf('val_predictions_%s', method), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## just to check validity of other layers, never used later on
    
    ## CAN-level bivariate scatterplots
    
    ## observed-observed
    if (method == params3$methods[1]) {  ## to run this only once as the observed values are only one
      x.axis.var <- results.shp.df@data$O_el_p95   ## x-axis is always the observed elev_p95 (main variable)
      for ( idx.targ in which(!unlist(params3$targ.names.sh) %in% "el_p95") ) {   ## loop over the rest of the resp. variables
        cmd <- sprintf("y.axis.var <- results.shp.df@data$O_%s", params3$targ.names.sh[idx.targ])   ## set as y-axis the other observed resp. variable
        eval(parse(text=cmd))
        fig.name.str <- file.path(BivCheck.CAN.subdir, sprintf("Bivariate_ObsObs_scatter_elevp95_vs_%s.pdf", params3$targ.names.sh[idx.targ]), sep='')
        pdf(fig.name.str)
          theme_set(theme_gray(base_size = 18))
          plot.to.print <- plot_colorByDensity(x.axis.var, y.axis.var, xlab="obs. elev_p95 [m]", ylab=sprintf("obs. %s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ]), main=sprintf("R=%.3f", cor(x.axis.var, y.axis.var)))  ## add simple correlation 
          print(plot.to.print)
        dev.off()
      }
    }  
    
    ## observed-predicted
    x.axis.var <- results.shp.df@data$O_el_p95
    for ( idx.targ in which(!unlist(params3$targ.names.sh) %in% "el_p95") ) {
      cmd <- sprintf("y.axis.var <- results.shp.df@data$P_%s", params3$targ.names.sh[idx.targ])  ## set as y-axis the other predicted resp. variable
      eval(parse(text=cmd))
      fig.name.str <- file.path(BivCheck.CAN.subdir, sprintf("Bivariate_ObsPred_scatter_elevp95_vs_%s_%s.pdf", params3$targ.names.sh[idx.targ], method), sep='')
      pdf(fig.name.str)
        theme_set(theme_gray(base_size = 18))
        plot.to.print <- plot_colorByDensity(x.axis.var, y.axis.var, xlab="obs. elev_p95 [m]", ylab=sprintf("pred. %s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ]), main=sprintf("R=%.3f", cor(x.axis.var, y.axis.var)))
        print(plot.to.print)
      dev.off()
    }
    
    ## predicted-predicted  
    x.axis.var <- results.shp.df@data$P_el_p95   ## x-axis is always the predicted elev_p95
    for ( idx.targ in which(!unlist(params3$targ.names.sh) %in% "el_p95") ) {
      cmd <- sprintf("y.axis.var <- results.shp.df@data$P_%s", params3$targ.names.sh[idx.targ])  ## set as y-axis the other predicted resp. variable
      eval(parse(text=cmd))
      fig.name.str <- file.path(BivCheck.CAN.subdir, sprintf("Bivariate_PredPred_scatter_elevp95_vs_%s_%s.pdf", params3$targ.names.sh[idx.targ], method), sep='')
      pdf(fig.name.str)
        theme_set(theme_gray(base_size = 18))
        plot.to.print <- plot_colorByDensity(x.axis.var, y.axis.var, xlab="pred. elev_p95 [m]", ylab=sprintf("pred. %s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ]), main=sprintf("R=%.3f", cor(x.axis.var, y.axis.var)))
        print(plot.to.print)
      dev.off()
    }
    
    ## residual-residual
    x.axis.var <- results.shp.df@data$R_el_p95   ## x-axis is always the residuals of elev_p95
    ## to have same axis limits... 
    if (method == params3$methods[1]) {  ## ...only at the 1st round of the method loop (for RF)...
      lims.table.x <- matrix(nrow = length(params3$targ.names.sh), ncol=2)  ## ...initialize limits for x and y axis (2-columns matrix) and then...
      lims.table.y <- lims.table.x
    }
    for ( idx.targ in which(!unlist(params3$targ.names.sh) %in% "el_p95") ) {
      cmd <- sprintf("y.axis.var <- results.shp.df@data$R_%s", params3$targ.names.sh[idx.targ])
      eval(parse(text=cmd))
      if (method == params3$methods[1]) {  ## ...compute only once, for RF, the range to be plotted and reuse the same for the YAI
        lims.table.x[idx.targ, ] <- c(quantile(x.axis.var, 0.0001, names=F), quantile(x.axis.var, 0.9999, names=F))
        lims.table.y[idx.targ, ] <- c(quantile(y.axis.var, 0.0001, names=F), quantile(y.axis.var, 0.9999, names=F))
      }
      fig.name.str <- file.path(BivCheck.CAN.subdir, sprintf("Bivariate_ResidResid_scatter_elevp95_vs_%s_%s.pdf", params3$targ.names.sh[idx.targ], method), sep='')
      pdf(fig.name.str)
        theme_set(theme_gray(base_size = 18))
        plot.to.print <- plot_colorByDensity(x.axis.var, y.axis.var, xlim=lims.table.x[idx.targ, ], ylim=lims.table.y[idx.targ, ], xlab="resid. elev_p95 [m]", ylab=sprintf("resid. %s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ]), main=sprintf("R=%.3f", cor(x.axis.var, y.axis.var)))
        if (params3$targ.names.sh[idx.targ] == "p_1r_2m") {   ## if the comparison is with cover_2m (p_1r_2m) add the vertical lines showing incoherent predictions (beyond params3$distort.thresh.<targ>)
          plot.to.print <- plot.to.print + geom_vline(xintercept=c(-params3$distort.thresh.elevp95, params3$distort.thresh.elevp95)) + geom_hline(yintercept=c(-params3$distort.thresh.cover2m, params3$distort.thresh.cover2m))
        }
        print(plot.to.print)
      dev.off()
    }

    ## Correlation and corr difference matrices for observed/predicted values on validation set
    
    if (method == params3$methods[1]) {  ## compute observed correlation matrix only once
      obs.data.matrix <- results.shp.df@data[ , seq(1,ncol(results.shp.df)-3,3)]  ## select every 3rd column until 3rd to last column to only use observed and predicted values
      colnames(obs.data.matrix) <- params3$targ.names.plots
      obs.corr.matrix <- cor(obs.data.matrix)
      write.csv(obs.corr.matrix, file = file.path(base_results_dir, "Val_obs_corr_matrix.csv", sep = ''))  ## save csv and tex with name telling it is the observed correlation matrix
      print( xtable(obs.corr.matrix, digits=rep(2, length(ncol(obs.corr.matrix)))), file = file.path(base_results_dir, "Val_obs_corr_matrix.tex", sep = '') )  
    }
    pred.data.matrix <- results.shp.df@data[ , seq(2,ncol(results.shp.df)-3,3)]  ## select every 3rd column starting from 2nd
    colnames(pred.data.matrix) <- params3$targ.names.plots
    pred.corr.matrix <- cor(pred.data.matrix)
    write.csv(pred.corr.matrix, file = file.path(base_results_dir, sprintf("Val_pred_%s_corr_matrix.csv", method), sep = ''))  ## save csv and tex with name telling it is the predicted correlation matrix for a given method
    cmd <- sprintf( "print( xtable(pred.corr.matrix, digits=rep(2, length(ncol(pred.corr.matrix)))), file = file.path(base_results_dir, \"Val_pred_%s_corr_matrix.tex\", sep = '') )", method)
    eval(parse(text=cmd))
    diff.corr.matrix <- pred.corr.matrix-obs.corr.matrix
    write.csv(diff.corr.matrix, file = file.path(base_results_dir, sprintf("Val_diff_%s_corr_matrix.csv", method), sep = ''))
    cmd <- sprintf( "print( xtable(diff.corr.matrix, digits=rep(2, length(ncol(diff.corr.matrix)))), file = file.path(base_results_dir, \"Val_diff_%s_corr_matrix.tex\", sep = '') )", method)
    eval(parse(text=cmd))
    

#### SUMMARY TABLES AND PLOTS BY CATEGORY ---------------------------------------------------------------
    
    ## round data to 1, 2 or 3 digits depending on the column
    CAN.stats.df.rounded <- round(CAN.stats.df, 2)
    CAN.stats.df.rounded[, c("rsq", "ac_uns", "ac_sys")] <- round(CAN.stats.df[, c("rsq", "ac_uns", "ac_sys")], 3)
    CAN.stats.df.rounded[, "rmsdpct"] <- round(CAN.stats.df[, "rmsdpct"], 1)
    
    ## add nr of samples per category
    ECO.sorting.idx <- match(colnames(ECO.stats.df), as.character(names(TRN.samples.per.Ecozone))) ## to match order of list of Ecozones
    ECO.stats.df <- rbind(round(ECO.stats.df, 3), as.integer(TRN.samples.per.Ecozone)[ECO.sorting.idx], VAL.samples.per.Ecozone[ECO.sorting.idx])
    rownames(ECO.stats.df)[length(params3$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
    UTM.sorting.idx <- match(as.character(zone.nrs$UTMzone), as.character(names(TRN.samples.per.UTMzone)))  ## to match order of list of UTM zones
    UTM.stats.df <- rbind(round(UTM.stats.df, 3), TRN.samples.per.UTMzone[UTM.sorting.idx], VAL.samples.per.UTMzone[UTM.sorting.idx])
    rownames(UTM.stats.df)[length(params3$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
    CNG.sorting.idx <- match(params3$Ch_attr.labels, as.character(names(TRN.samples.per.Ch_attr)))
    CNG.stats.df <- rbind(round(CNG.stats.df, 3), TRN.samples.per.Ch_attr[CNG.sorting.idx], VAL.samples.per.Ch_attr[CNG.sorting.idx])
    rownames(CNG.stats.df)[length(params3$targ.names.lg)+1:2] <- c("Nr. TRN samples", "Nr. VAL samples")
    
    ## stack results in df stats3
    cmd <- sprintf('stats3$%s$CAN.stats.df <- CAN.stats.df.rounded', method)
    eval(parse(text=cmd))
    cmd <- sprintf('stats3$%s$ECO.stats.df <- ECO.stats.df', method)
    eval(parse(text=cmd))
    cmd <- sprintf('stats3$%s$UTM.stats.df <- UTM.stats.df', method)
    eval(parse(text=cmd))
    cmd <- sprintf('stats3$%s$CNG.stats.df <- CNG.stats.df', method)
    eval(parse(text=cmd))
    
    ## save csv and tex files
    cmd <- sprintf( "write.csv( stats3$%s$CAN.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_CAN_stats.csv\", sep = '') )", method, method )
    eval(parse(text=cmd))
    cmd <- sprintf( "print( xtable(stats3$%s$CAN.stats.df, digits=c(0,2,2,2,2,2,2,2,2,3,2,1,3,3,2)), file = file.path(base_results_dir, \"RES_MA_%s_CAN_stats.tex\", sep = '') )", method, method )
    eval(parse(text=cmd))
    
    cmd <- sprintf( "write.csv( stats3$%s$ECO.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_ECO_%s.csv\", sep = '') )", method, method, params3$metrics.MS )
    eval(parse(text=cmd))
    cmd <- sprintf( "print( xtable(stats3$%s$ECO.stats.df, digits=c(0,rep(3, length(params3$sampled.ecozones)))), file = file.path(base_results_dir, \"RES_MA_%s_ECO_%s.tex\", sep = '') )", method, method, params3$metrics.MS)
    eval(parse(text=cmd))
    
    cmd <- sprintf( "write.csv( stats3$%s$UTM.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_UTM_%s.csv\", sep = '') )", method, method, params3$metrics.MS )
    eval(parse(text=cmd))
    cmd <- sprintf( "print( xtable(stats3$%s$UTM.stats.df, digits=c(0,rep(3, length(paramsGL$zones)))), file = file.path(base_results_dir, \"RES_MA_%s_UTM_%s.tex\", sep = '') )", method, method, params3$metrics.MS)
    eval(parse(text=cmd))
    
    cmd <- sprintf( "write.csv( stats3$%s$CNG.stats.df, file = file.path(base_results_dir, \"RES_MA_%s_CNG_%s.csv\", sep = '') )", method, method, params3$metrics.MS )
    eval(parse(text=cmd))
    cmd <- sprintf( "print( xtable(stats3$%s$CNG.stats.df, digits=c(0,rep(3, 4))), file = file.path(base_results_dir, \"RES_MA_%s_CNG_%s.tex\", sep = '') )", method, method, params3$metrics.MS)
    eval(parse(text=cmd))
    
    ## Rsq barplots by category
    
    targs.for.barplots <- c("elev_p95", "cover_2m", "ag_biomass")  ## only for 3 key resp. variables
    
    ## by Ecozone
    data.to.plot <- t(ECO.stats.df[targs.for.barplots, ])  ## select the 3 rows and transpose them 
    data.melted <- melt(data.to.plot)   ## melt them to accomodate ggplot needs for barplots
    colnames(data.melted) <- c("Ecozone", "Response_var", "R2")
    fig.name.str <- file.path(Assess.CAN.subdir, sprintf("Barplot_ECO_%s_%s.pdf", method, params3$metrics.MS), sep='')
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
                       legend.position="top",
                       legend.direction="horizontal") +
                 scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8)) +
                 ylab("R^2")
    print(barplot)
    dev.off()
    
    
    ## by UTM zone
    data.to.plot <- t(UTM.stats.df[targs.for.barplots, ])
    data.melted <- melt(data.to.plot)
    colnames(data.melted) <- c("UTMzone", "Response_var", "R2")
    fig.name.str <- file.path(Assess.CAN.subdir, sprintf("Barplot_UTM_%s_%s.pdf", method, params3$metrics.MS), sep='')
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
    fig.name.str <- file.path(Assess.CAN.subdir, sprintf("Barplot_CNG_%s_%s.pdf", method, params3$metrics.MS), sep='')
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
  
  }  ## end for method in params3$methods
  
  stats.file = file.path(base_results_dir, 'stats3.Rdata', fsep = .Platform$file.sep) 
  save(stats3, file = stats.file)

}  ## end if params3$run.SG

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Prog3, total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)
