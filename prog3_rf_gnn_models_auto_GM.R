########################################################
## Saskatchewan random forest and gnn models
## response variables from lidar plots
## explanatory variables from landsat spectral indices landsat change metrics and topographic indices
## explanatory variables extracted for lidar plots
## RandomForest analysis and variable selection
## GNN CCA analysis and variable selection 
## save RF GNN models and training data as an Rdata file
#######################################################
##----------------------
## TO DO
##----------------------

## STILL TO DO
# - after selection of the final set of variables save final names, etc.

## SOLVED
# -V write script between prog2 and prog3 to group data by ecozone rather than by UTM zone -- no, we run the analysis at the national level with one single model (long, lat & other trends should account)
# -V change loop over ecozones instead of over UTM zones -- no need to do it as we use one single model valid for all Can territory


##----------------------
## READS
##----------------------
# -

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

##----------------------------
## SCRIPT SPECIFIC PARAMETERS
##----------------------------

## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end
nr.clusters = length(zones)

##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("rgdal",
                      "raster",
                      "sp",
                      "randomForest",
                      "rgl",
                      "yaImpute",
                      "vegan",
                      "snow",
                      "lubridate", 
                      "doParallel", 
                      "foreach"
                      # "profvis",
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  cmd=sprintf('library(%s)', pack)
  eval(parse(text=cmd))
}

tic <- proc.time() # start clocking global time


print('Prog3, modeling') 
wkg_dir = base_wkg_dir
results_dir = base_results_dir

# setwd(wkg_dir)

# ## set basepath
# basepath<-"F:/SK_nn_forest/GeordieWorking/wkg/new_model"
# setwd(basepath)


##load plot training and validation data for lidar metrics and derived structural attributes (old "str", for structural)
##subset to only training data
Y.trn.val <- read.csv(file.path(base_wkg_dir, "lidar_metrics_mean_training_validation.csv", fsep = .Platform$file.sep))
Y.trn.val <- Y.trn.val[,-1] # remove 1st column X

##load explanatory variables extracted for training plots (old "env" for environment)
X.trn.val <- read.csv(file.path(base_wkg_dir, "poly_training_validation_exvars_extract.csv", fsep = .Platform$file.sep))

##change POLY250ID to FCID in str data
names(Y.trn.val)[names(Y.trn.val)=="POLY250ID"] <- "FCID"
names(Y.trn.val)

##remove unneeded vars in env (X.trn.val) data
##remove years since disturbance retaiing only calendar year of firest greatest and last disturbance and add in the end the TV variable for selection
names(X.trn.val)
X.trn.val <- X.trn.val[,c(2,6:10,12:17,23:25,27,29,32:34, 4)]

##extract core str response variables using in random forest (RF) analysis
##core str response variables are elev_mean elev_stddev elev_p95 elev_cv percentage_first_returns_above_2m percentage_first_returns_above_mean
##imputation map only contains single band of plot FCIDs
##extract key variables that get attached to map after imputation
Y.trn.val.sub <- Y.trn.val[,c("FCID","elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean", "loreys_height", "basal_area", "gross_stem_volume", "total_biomass", "TV")]
Y.trn.val.attach <- Y.trn.val[,c("FCID","elev_mean","elev_stddev","elev_p95","elev_cv","percentage_first_returns_above_2m","percentage_first_returns_above_mean",
                                    "mean_tree_height","dominant_tree_height","loreys_height","basal_area","gross_stem_volume","foliage_biomass",
                                    "branch_biomass","crown_biomass","bark_biomass","wood_biomass","stem_biomass","total_biomass", "TV")]

names(Y.trn.val)
names(Y.trn.val.sub)

## function to convert row names of dataframes to FCID
RownamesToFCID <- function(indataframe){
  rownames(indataframe) <- indataframe[,"FCID"]
  indataframe <- indataframe[,2:ncol(indataframe)]
  return(indataframe)
}
X.trn.val <- RownamesToFCID(X.trn.val)
Y.trn.val.sub <- RownamesToFCID(Y.trn.val.sub)
Y.trn.val.attach <- RownamesToFCID(Y.trn.val.attach)

X.trn <- subset(X.trn.val, TV == "TRAINING")
Y.trn.sub <- subset(Y.trn.val.sub, TV == "TRAINING")
Y.trn.attach <- subset(Y.trn.val.attach, TV == "TRAINING")
#   X.val <- subset(X.trn.val, TV == "VALIDATION")
#   Y.val.sub <- subset(Y.trn.val.sub, TV == "VALIDATION")
#   Y.val.attach <- subset(Y.trn.val.attach, TV == "VALIDATION")

# remove column with TV label
X.trn <- X.trn[,-ncol(X.trn)]
Y.trn.sub <- Y.trn.sub[,-ncol(Y.trn.sub)]
Y.trn.attach <- Y.trn.attach[,-ncol(Y.trn.attach)]
#   X.val <- X.val[,-ncol(X.val)]
#   Y.val.sub <- Y.val.sub[,-ncol(Y.val.sub)]
#   Y.val.attach <- Y.val.attach[,-ncol(Y.val.attach)]

##check for rows of data with values totaling zero WHY???
min(apply(X.trn, MARGIN = 1, sum)) ## all good
min(apply(Y.trn.sub, MARGIN = 1, sum)) ## all good
min(apply(Y.trn.attach, MARGIN = 1, sum)) ## all good

## sort dataframes by rownames
X.trn <- X.trn[order(row.names(X.trn)),]
Y.trn.sub <- Y.trn.sub[order(row.names(Y.trn.sub)),]
Y.trn.attach <- Y.trn.attach[order(row.names(Y.trn.attach)),]

##final dim check and names check 
dim(X.trn)
dim(Y.trn.sub)
dim(Y.trn.attach)
names(X.trn)
names(Y.trn.sub)
names(Y.trn.attach)

########################################################################
## Random Forest Model with 500 trees
set.seed(2010)
names(X.trn)

#   yai.rf.500 <- yai(x=X.trn[,c("TCB","TCG","TCW","TCA","TCD",
#                                "PreChange_persistence","PreChange_magnitude","PreChange_evolution",
#                                "PostChange_persistence","PostChange_magnitude","PostChange_evolution",
#                                "Years_Since_Greatest_Change","Change_rate", "Change_persistence","Change_magnitude",
#                                 "Elevation","Slope","TWI","TSRI")],y=str, method="randomForest", k = 1, ntree=500)

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

#   yai.rf.250 <- yai(x=X.trn.val[,c("TCB","TCG","TCW","TCA","TCD",
#                                    "PreChange_persistence","PreChange_magnitude","PreChange_evolution",
#                                    "PostChange_persistence","PostChange_magnitude","PostChange_evolution",
#                                    "Years_Since_Greatest_Change","Change_rate", "Change_persistence","Change_magnitude",
#                                    "Elevation","Slope","TWI","TSRI")],y=str, method="randomForest", k = 1, ntree=250)
#   

yai.rf.250 <- yai(x=X.trn, y=Y.trn.sub, method="randomForest", k = 1, ntree=250)


##variable importance plot
yaiVarImp.250 <- yaiVarImp(yai.rf.250, nTop =0, plot=TRUE)
str(yaiVarImp.250)
par(las=1)
par(mar=c(4,15,1,1))
#   pdf("c:/deleteme/test3.pdf")
boxplot(yaiVarImp.250, horizontal = TRUE, xlab = "Scaled Variable Importance")

#########################################################################
#   ## GNN CCA model 
#   yai.gnn.sk <- yai(x=X.trn.val[,c("TCB","TCG","TCW","TCA","TCD",
#                                 "PreChange_persistence","PreChange_magnitude","PreChange_evolution",
#                                 "PostChange_persistence","PostChange_magnitude","PostChange_evolution",
#                                 "Years_Since_Greatest_Change","Change_rate", "Change_persistence","Change_magnitude",
#                                 "Elevation","Slope","TWI","TSRI")],y=str, method="gnn", k = 1)
#   


##########################################################################
##save data and models as an RData file for future easy loading
## save rf and cca models
## save plot data both 6 response vars used in models and anncilary plot vars to be attached to maps
## save predcitor vars associated with plots
save(yai.rf.500, yai.rf.250, X.trn, Y.trn.sub, Y.trn.attach, file = "yai_rf_gnn_sk.RData")

# clock UTM zone time
temp.toc <- proc.time()-temp.tic[3]
print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))



# stopCluster(cl)

# clock global time
toc <- proc.time()-tic[3]
print(paste("Prog3, total elapsed time:",seconds_to_period(toc[3])))