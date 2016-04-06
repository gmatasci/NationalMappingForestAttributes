########################################################
## response variables from lidar plots
## explanatory variables from landsat spectral indices landsat change metrics and topographic indices
## explanatory variables extracted for lidar plots
## RandomForest analysis and variable selection
#######################################################
##----------------------
## TO DO
##----------------------

## STILL TO DO
# - check plots for change vars with new averages and trends
# - check plot section with new structure 
# - ANOVA style qplot for catgorical vars
# - ggpair each type vs the targets: fix loop for dynamic plots (consider copying the same 4 times)
# - implement model selection based on Impute_Vs_predict.R
# - build model with proprer data
# - separate script (not to be run everytime) for model selection called "prog3a_" to see if there
#   are improvements in predicting response variables on at a time and not as a single multivariate Y: 
#   this means comparing yaimpute rf, yaimpute gnn, rf, svm, etc.
# - yaImpute uses a RF with a nrTrees shared across all Ys, so actual nrTrees is nrTrees/nrYs
# - ############# = TODEL

## SOLVED
# -V write script between prog2 and prog3 to group data by ecozone rather than by UTM zone -- no, we run the analysis at the national level with one single model (long, lat & other trends should account)
# -V change loop over ecozones instead of over UTM zones -- no need to do it as we use one single model valid for all Can territory
# -V multi.hist() per var by UTM zone -- no, only works for dataframes with same nr of columns (UTM zones have different nr of samples)


##----------------------
## READS
##----------------------
# - "lidar_metrics_mean_training_validation.csv": (from prog1) csv table with all observed LiDAR and forest attributes (Y) for selected TRN and VAL samples (3x3 polygons with average plot values or just single plot value) for all UTM zones
# - "poly_training_validation_exvars_extract.csv":  (from prog2) csv table with explanatory variables (X) for selected TRN and VAL samples (3x3 plots or single plot polygons with weighted average pixels values)

##----------------------
## WRITES
##----------------------
# - "yai_rf_gnn_sk.RData": R data file with yaimpute models and training set to use later for mapping (prediction/imputation) on the grid

#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

rm(list=ls()) # clear all variables

##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

source("Functions_NatImp.R")


##----------------------------
## SCRIPT SPECIFIC PARAMETERS
##----------------------------
params3a <- list()

params3a$targ.names.lg <- list("elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean", "loreys_height", "basal_area", "gross_stem_volume", "total_biomass")
params3a$targ.names.sh <- list("el_m","el_std","el_p95","el_cv","pct_1r_ab2","pct_1r_abm", "loreys_h", "basal_a", "stem_v", "tot_biom")

params3a$targ.key.names.lg <- list("percentage_first_returns_above_2m", "total_biomass")
params3a$targ.key.names.sh <- list("pct_1r_ab2", "tot_biom")

params3a$targ.attach.names.lg <- list("elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean",
                                "mean_tree_height","dominant_tree_height","loreys_height","basal_area","gross_stem_volume","foliage_biomass",
                                "branch_biomass","crown_biomass","bark_biomass","wood_biomass","stem_biomass","total_biomass")


# params3a$subs.factor <- 100
params3a$subs.factor <- 20

params3a$groups.to.plot <- list('Bands', 'TCcomps', 'Change', 'PrePostChange', 'GreatestLastFirstChange', 'ChangeAttribution', 'Trends')

param_file_prog3a = file.path(base_wkg_dir, 'AllUTMzones_params3a.Rdata', fsep = .Platform$file.sep) 
save(params3a, file = param_file_prog3a)


##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("ggplot2", 
                      "GGally",
                      "psych",
                      "randomForest",
                      "rgl",
                      "yaImpute",
                      "vegan",
                      "snow",
                      "lubridate", 
                      "doParallel", 
                      "foreach"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

tic <- proc.time() # start clocking global time


print('Prog3, modeling') 

##load plot training and validation data for lidar metrics and derived structural attributes (old "str", for structural)
##subset to only training data
Y.trn.val <- read.csv(file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", fsep = .Platform$file.sep))
# Y.trn.val$POLY250ID <- NULL

##load explanatory variables extracted for training plots (old "env" for environment)
X.trn.val <- read.csv(file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", fsep = .Platform$file.sep))

##change POLY250ID to FCID in Y data
names(Y.trn.val)[names(Y.trn.val)=="POLY250ID"] <- "FCID"

## check for order of sample points in the two datasets X and Y
if ( !identical(X.trn.val$FCID, Y.trn.val$FCID) ) {
  stop("Prog3a: order of FCIDs is different in lidar_metrics_mean_training_validation.csv wrt poly_training_validation_exvars_extract.csv")
}

##------------------------
## DESCRIPTIVE PLOTS
##------------------------

## create indices for subsampling 
plot.idx <- sample(1:nrow(X.trn.val), round(nrow(X.trn.val)/params3a$subs.factor))

## Canadian level (merged UTM zones) bivariate scatterplots

Y.data.to.plot <- Y.trn.val[plot.idx, unlist(params3a$targ.key.names.lg)] ## select Y columns to plot
colnames(Y.data.to.plot) <- params3a$targ.key.names.sh

gr.idx <- 1
for (pred.gr.names in params3a$groups.to.plot) {
  gr.name <- unlist(params3a$groups.to.plot[gr.idx])
  
  if (gr.name == 'ChangeAttribution') {
    str <- file.path(base_figures_dir, sprintf("CAN_boxplots_keyTarg_VS_%s.pdf", gr.name), sep='')
    pdf(str)
      par(mfrow=c(1,2))
      describeBy(Y.data.to.plot$pct_1r_ab2, group = X.trn.val[, "Ch_attr"])
      qplot(Ch_attr, Y.data.to.plot$pct_1r_ab2, data=X.trn.val, geom="boxplot")
      describeBy(Y.data.to.plot$tot_biom, group = X.trn.val[, "Ch_attr"])
      qplot(Ch_attr, Y.data.to.plot$tot_biom, data=X.trn.val, geom="boxplot")
      title(sprintf("%s VS pct_1ret_ab_2m and tot_biom", gr.name), outer=TRUE)
    dev.off()
  } else if (gr.name == 'Bands') {
    pred.names <- unlist(params2$pred.names.sh$bands)
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  } else if (gr.name == 'TCcomps') {
    pred.names <- unlist(params2$pred.names.sh$TC)
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  } else if (gr.name == 'Change') {
    pred.names <- c("Ch_pers", "Ch_mag", "Ch_er")
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  } else if (gr.name == 'PrePostChange') {
    pred.names <- c("PreCh_pers", "PreCh_mag", "PreCh_er", "PostCh_pers", "PostCh_mag", "PostCh_er")
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  } else if (gr.name == 'GreatestLastFirstChange') {
    pred.names <- c("GreatCh_yr", "FirstCh_yr", "LastCh_yr", "FirstCh_pers", "LastCh_pers")
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  } else if (gr.name == 'Trends') {
    pred.names <- unlist(params2$pred.names.sh$trends)
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  }
  
  X.data.to.plot <- X.trn.val[plot.idx, pred.names]
  data.to.plot <- cbind(X.data.to.plot, Y.data.to.plot)
  str <- file.path(base_figures_dir, sprintf("CAN_pairsPanel_keyTarg_VS_%s.pdf", gr.name), sep='')
  pdf(str)
    pairs.panels(data.to.plot, pch=1, cex=cex.val, 
                 scale=T, density=F, ellipses=F, lm=T, jiggle=F, rug=F,
                 breaks=15, cex.labels=cex.labels.val, cex.cor=cex.cor.val,
                 main=sprintf("%s VS pct_1ret_ab_2m and tot_biom", gr.name))
  dev.off()
  
  gr.idx <- gr.idx+1
}


## Variable histograms by UTM zone

## predictors
X.data.to.plot <- X.trn.val[plot.idx, ]
for (pred in params2$var.names.short) {
  str <- file.path(base_figures_dir, sprintf("UTMzonesHist_%s.pdf", pred), sep='')
  pdf(str)
  par(mfrow=c(5,4))
  for (z in 1:length(paramsGL$zones)) {
    zone <- paramsGL$zones[z]
    zone.nr <- substr(zone, 4, nchar(zone))
    zone.idx <- which(X.trn.val$UTMzone == zone.nr)
    hist(X.trn.val[zone.idx, pred], main=sprintf("%s", zone))  # , breaks=20
  }
  dev.off()
}
    
## targets
Y.data.to.plot <- Y.trn.val[plot.idx, ]
for (targ in params3a$targets) {
  str <- file.path(base_figures_dir, sprintf("UTMzonesHist_%s.pdf", targ), sep='')
  pdf(str)
  par(mfrow=c(5,4))
  for (z in 1:length(paramsGL$zones)) {
    zone <- paramsGL$zones[z]
    zone.nr <- substr(zone, 4, nchar(zone))
    zone.idx <- which(X.trn.val$UTMzone == zone.nr)
    hist(Y.data.to.plot[zone.idx, targ], main=sprintf("%s", zone))  # , breaks=20
  }
  dev.off()
}



##------------------------
## MODEL SELECTION
##------------------------

############# = TODEL
##remove unneeded vars in env (X.trn.val) data
##remove years since disturbance retaiing only calendar year of firest greatest and last disturbance and add in the end the TV variable for selection
# names(X.trn.val)
# X.trn.val <- X.trn.val[,c(2,6:10,12:17,23:25,27,29,32:34, 4)]
############# = TODEL

##remove unneeded vars in X.trn.val data
X.trn.val <- X.trn.val[, c('FCID', unlist(params2$pred.names.sh), 'TV')]

##extract core str response variables using in random forest (RF) analysis
##core str response variables are elev_mean elev_stddev elev_p95 elev_cv percentage_first_returns_above_2m percentage_first_returns_above_mean
##imputation map only contains single band of plot FCIDs
##extract key variables that get attached to map after imputation
Y.trn.val.sub <- Y.trn.val[, c("FCID", unlist(params3a$targ.names.lg), "TV")]
Y.trn.val.attach <- Y.trn.val[, c("FCID", unlist(params3a$targ.attach.names.lg), "TV")]

X.trn.val <- RownamesToFCID(X.trn.val)
Y.trn.val.sub <- RownamesToFCID(Y.trn.val.sub)
Y.trn.val.attach <- RownamesToFCID(Y.trn.val.attach)

X.trn <- subset(X.trn.val, TV == "TRAINING")
Y.trn.sub <- subset(Y.trn.val.sub, TV == "TRAINING")
Y.trn.attach <- subset(Y.trn.val.attach, TV == "TRAINING")
X.val <- subset(X.trn.val, TV == "VALIDATION")
Y.val.sub <- subset(Y.trn.val.sub, TV == "VALIDATION")
Y.val.attach <- subset(Y.trn.val.attach, TV == "VALIDATION")

# remove column with TV label
X.trn$TV <- NULL
Y.trn.sub$TV <- NULL
Y.trn.attach$TV <- NULL
X.val$TV <- NULL
Y.val.sub$TV <- NULL
Y.val.attach$TV <- NULL
  
############# = TODEL
##check for rows of data with values totaling zero WHY???
# min(apply(X.trn, MARGIN = 1, sum)) ## all good
# min(apply(Y.trn.sub, MARGIN = 1, sum)) ## all good
# min(apply(Y.trn.attach, MARGIN = 1, sum)) ## all good

#!!!!! TO CHECK WHY????? AND IN CASE SORT ALSO VAL
## sort dataframes by rownames
X.trn <- X.trn[order(row.names(X.trn)),]
Y.trn.sub <- Y.trn.sub[order(row.names(Y.trn.sub)),]
Y.trn.attach <- Y.trn.attach[order(row.names(Y.trn.attach)),]

##final dim check and names check 
# dim(X.trn)
# dim(Y.trn.sub)
# dim(Y.trn.attach)
# names(X.trn)
# names(Y.trn.sub)
# names(Y.trn.attach)
############# = TODEL

########################################################################
## Random Forest Model with 500 trees
set.seed(2010)

yai.rf.500 <- yai(x=X.trn, y=Y.trn.sub, method="randomForest", k = 1, ntree=500)

##variable importance plot
yaiVarImp.500 <- yaiVarImp(yai.rf.500, nTop =0, plot=TRUE)
str(yaiVarImp.500)
par(las=1)
par(mar=c(4,15,1,1))
  boxplot(yaiVarImp.500, horizontal = TRUE, xlab = "Scaled Variable Importance")
dev.off()

########################################################################
## Random Forest Model with 250 trees 
set.seed(2010)
names(X.trn.val)

yai.rf.250 <- yai(x=X.trn, y=Y.trn.sub, method="randomForest", k = 1, ntree=250)


##variable importance plot
yaiVarImp.250 <- yaiVarImp(yai.rf.250, nTop =0, plot=TRUE)
str(yaiVarImp.250)
par(las=1)
par(mar=c(4,15,1,1))
#   pdf("c:/deleteme/test3.pdf")
boxplot(yaiVarImp.250, horizontal = TRUE, xlab = "Scaled Variable Importance")

##########################################################################
##save data and models as an RData file for future easy loading
## save rf and cca models
## save plot data both 6 response vars used in models and anncilary plot vars to be attached to maps
## save predcitor vars associated with plots
save(yai.rf.500, yai.rf.250, X.trn, Y.trn.sub, Y.trn.attach, file = "yai_rf_gnn_sk.RData")


##------------------------
## MODEL ASSESSMENT
##------------------------





# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog3, total elapsed time:",seconds_to_period(toc[3])))