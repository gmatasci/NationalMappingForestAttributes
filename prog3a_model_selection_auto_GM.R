########################################################
## response variables from lidar plots
## explanatory variables from landsat spectral indices landsat change metrics and topographic indices
## explanatory variables extracted for lidar plots
## RandomForest analysis and variable selection
#######################################################
##----------------------
## TO DO
##----------------------

#### STILL TO DO ####
# - consider setting mtry=sqrt(nr var)
# - 10 run and report mean+-std?
# - subsample only for pair.plot and not for hist or rest
# - implement model selection based on Impute_Vs_predict.R
# - build model with proprer data
# - run with final data: comparison between predicting response variables one at a time and not as a single multivariate Y: 
#   this means comparing yaimpute rf, rf, etc.
# - yaImpute uses a RF with a nrTrees shared across all Ys, so actual nrTrees is nrTrees/nrYs
# - check plots for change vars with new averages and trends
# - make sure Ch_attr is interpreted as a categorical/factor variable by RF and that all the factor values are found in TRN (no new values in VAL), otherwise errors (http://stats.stackexchange.com/questions/29446/random-forest-and-new-factor-levels-in-test-set)
# - ############# = TODEL
# - check to see sampled values based on "FCID" or "POLY250ID" from X or Y dataframes loaded here linked to field "POLY250ID" of <UTMzone>_cpt_poly_250m_training_validation.shp (saved in /BAP_Imputation_working/wkg/<UTMzone>/)
#   then related to <UTMzone>_plot_inventory_attributes2.csv by "unique_id"

#### SOLVED ####
# -V write script between prog2 and prog3 to group data by ecozone rather than by UTM zone -- no, we run the analysis at the national level with one single model (long, lat & other trends should account)
# -V change loop over ecozones instead of over UTM zones -- no need to do it as we use one single model valid for all Can territory
# -V multi.hist() per var by UTM zone -- no, only works for dataframes with same nr of columns (UTM zones have different nr of samples)
# -V ANOVA style qplot for categorical vars -- boxplot of key targets by category for CAN and barplot by UTM zone
# -V ggpair each type vs the targets: fix loop for dynamic plots (consider copying the same 4 times) -- done with pairs.panels()
# -V check plots section with new structure
# -V test density=T for pair.panels -- no, this just adds a density curve to the histogram in the diagonal
# -V global correlation matrix without subsampling to see redundant variables -- added with filters for high and low correlations




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

print('Prog3a: descriptive plotting and model selection') 

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

params3a$min.year.of <- 1984

params3a$subs.factor <- 60

params3a$Ch_attr.classes <- c("No change", "Fire", "Harvesting", "Non-stand replacing", "Road", "Unclassified")

params3a$groups.to.plot <- list('Bands', 'TCcomps_VI', 'ChangeAttribution', 'Years_Topo_Trends')

param_file_prog3a = file.path(base_wkg_dir, 'AllUTMzones_params3a.Rdata', fsep = .Platform$file.sep) 
save(params3a, file = param_file_prog3a)


##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("ggplot2",
                      "gridExtra",
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
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

tic <- proc.time() ## start clocking global time

## load plot training and validation data for lidar metrics and derived structural attributes (old "str", for structural)
Y.trn.val <- read.csv(file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", fsep = .Platform$file.sep))

## load explanatory variables extracted for training plots (old "env" for environment)
X.trn.val <- read.csv(file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", fsep = .Platform$file.sep))

## check for order of sample points in the two datasets X and Y
if ( !identical(X.trn.val$FCID, Y.trn.val$FCID) ) {
  stop("Prog3a: order of FCIDs is different in lidar_metrics_mean_training_validation.csv wrt poly_training_validation_exvars_extract.csv")
}

##---------------------------
## DESCRIPTIVE STATS & PLOTS
##---------------------------

CAN.subdir <- file.path(base_figures_dir, "Relations_CAN_level", fsep = .Platform$file.sep)
UTMzone.subdir <- file.path(base_figures_dir, "Histograms_UTMzone_level", fsep = .Platform$file.sep)

if (! file.exists(CAN.subdir)){dir.create(CAN.subdir, showWarnings = F, recursive = T)}
if (! file.exists(UTMzone.subdir)){dir.create(UTMzone.subdir, showWarnings = F, recursive = T)}

## global correlation matrix to see redundant variables
XXXX
XXXX Split to highlight high corr among predictors and low corr with response
XXXX

params3a$corr.high <- 0.95
params3a$corr.low <- 0.2
Y.trn.val.corr <- Y.trn.val[, unlist(params3a$targ.names.lg)]
colnames(Y.trn.val.corr) <- params3a$targ.names.sh
cont.idx <- !unlist(params2$pred.names.sh) %in% c(params2$pred.names.sh$yrs, params2$pred.names.sh$cngattr) ## subset to continuous variables
data.matrix <- cbind(X.trn.val[, unlist(params2$pred.names.sh)[cont.idx]], Y.trn.val.corr)
corr.matrix <- cor(data.matrix)
corr.matrix.high.low <- corr.matrix
corr.matrix.high.low[abs(corr.matrix) < params3a$corr.high & abs(corr.matrix) > params3a$corr.low] <- NA



## create indices for subsampling 
set.seed(2010)

plot.idx.init <- sample(1:nrow(X.trn.val), round(nrow(X.trn.val)/params3a$subs.factor))
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
    X.data.to.plot <- as.data.frame(X.trn.val[plot.idx, pred.names])
    colnames(X.data.to.plot) <- pred.names
    X.data.to.plot$Ch_attr <- as.factor(X.data.to.plot$Ch_attr)
    df.to.plot <- cbind(X.data.to.plot, Y.data.to.plot)
    df.to.plot$Ch_attr <- as.factor(df.to.plot$Ch_attr)
    
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
  } else if (gr.name == 'Years_Topo_Trends') {
    pred.names <- c(unlist(params2$pred.names.sh$yrs), unlist(params2$pred.names.sh$topo), unlist(params2$pred.names.sh$trends))
    cex.val <- 1
    cex.labels.val <- 1
    cex.cor.val <- 0.8
  }
  
  ## Canadian level (merged UTM zones) bivariate scatterplots with linear model R on it
  X.data.to.plot <- X.trn.val[plot.idx, pred.names]
  data.to.plot <- cbind(X.data.to.plot, Y.data.to.plot)
  if (gr.name == 'Years_Topo_Trends') {
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


## Variable histograms by UTM zone

## continuous var
cont.var.idx <- !(colnames(X.trn.val) %in% "Ch_attr")
data.to.plot <- cbind(X.trn.val[plot.idx, cont.var.idx], Y.data.to.plot)
data.to.plot$GreatCh_yr[data.to.plot$GreatCh_yr == 0] <- params3a$min.year.of

# XXXXXXXXXXXXX
# XXXXXXXXXXXXX
# - below change 6 bc columns have changed
# XXXXXXXXXXXXX
# XXXXXXXXXXXXX

for (v in 6:ncol(data.to.plot)) {
  var <- colnames(data.to.plot)[v]
  nr.of.bins <- 12
  my.hist.lims <- c(min(data.to.plot[, v]), quantile(data.to.plot[, v], 0.99, names=F))
  my.breaks <- c(seq(min(data.to.plot[, v]), quantile(data.to.plot[, v], 0.99, names=F), l=nr.of.bins+1), max(data.to.plot[, v])+1)
  str <- file.path(UTMzone.subdir, sprintf("UTMzonesHist_%s.pdf", var), sep='')
  pdf(str)
  par(mfrow=c(4,5), oma = c(5,4,2,2) + 0.1, mar = c(2,2,2,2) + 0.1)
  for (z in 1:length(paramsGL$zones)) {
    zone <- paramsGL$zones[z]
    zone.nr <- substr(zone, 4, nchar(zone))
    zone.idx <- which(data.to.plot$UTMzone == zone.nr)
    hist(data.to.plot[zone.idx, v], main=sprintf("%s", zone), freq=T, breaks=my.breaks, xlim=my.hist.lims, xlab=NULL, ylab=NULL, col='lightgreen')  # , breaks=20
  }
  title(var, outer=T)
  dev.off()
}
    
## categorical vars
data.to.plot <- cbind(X.trn.val[plot.idx, 1:5], X.trn.val[plot.idx, "Ch_attr"])
var <- "Ch_attr"
str <- file.path(UTMzone.subdir, sprintf("UTMzonesHist_%s.pdf", var), sep='')
pdf(str)
par(mfrow=c(4,5), oma = c(5,4,2,2) + 0.1, mar = c(0,0,2,2) + 0.1)
for (z in 1:length(paramsGL$zones)) {
  zone <- paramsGL$zones[z]
  zone.nr <- substr(zone, 4, nchar(zone))
  zone.idx <- which(data.to.plot$UTMzone == zone.nr)
  barplot(table(data.to.plot[, 6]), col='lightgreen')
}
title("Ch_attr", outer=T)
dev.off()













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

####      
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
####        



# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog3, total elapsed time:",seconds_to_period(toc[3])))