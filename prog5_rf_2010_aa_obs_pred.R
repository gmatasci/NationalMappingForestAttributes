################################################################
## saskatchewan imputation mapping 
## accuracy assessment of 2010 imputation map on validaton plots
## accuracy assessments based on observed versus predicted values
## does not include assessments in relation to magnitude of disturbance or time since disturbance
##----------------------
## TO DO
##----------------------
# - 

##----------------------
## READS
##----------------------
# - "sk_rf_2010_fcid.img": imputed IDs over all grid (30 m resolution) 
# - "sk_str_attach.csv": training Ys to attach to the predicted (imputed) IDs to get maps
# - "lidar_metrics_mean_training_validation.csv": csv table with all observed LiDAR and forest attributes (Ys) for selected TRN and VAL samples (3x3 polygons with average plot values)
# - "<UTMzone>_poly_250m_training_validation": polygon shp with training and validation polygons, subsetted to keep only validation polygons and used to weighted-average predicted (attached data to imputed IDs) LiDAR and forest attributes


##----------------------
## WRITES
##----------------------
# - "rf_2010_validation_obs_pred.csv": observed (coming from "lidar_metrics_mean_training_validation.csv") and 
#                                      predicted (coming from the merge of mapped FCIDs "sk_rf_2010_fcid.img" and "sk_str_attach.csv", then weighted-averaged over validation polygons) 
#                                      values on validation polygons 
# - "rf_2010_aa_metrics.csv": table with regression metrics assessing the predictions for all Ys 


#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

rm(list=ls()) # clear all variables

##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.R"
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
                      "ggplot2",
                      "grid",
                      "lubridate", 
                      "doParallel", 
                      "foreach"
                      # "profvis",
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

tic <- proc.time() # start clocking global time

# cl <- makeCluster(nr.clusters)
# registerDoParallel(cl)
# getDoParWorkers()
# foreach (z = 1:length(zones), .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
for (z in 1:length(zones)) {
  
  temp.tic <- proc.time() # start clocking time for each UTM zone
  
  zone = zones[z]
  print(paste('Validation on' ,zone))   # converts its arguments (via as.character) to character strings, and concatenates them (separating them by the string given by sep or by a space by default)
  wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)
  results_dir = file.path(base_results_dir, zone, fsep = .Platform$file.sep)
  
  setwd(wkg_dir)
  
  ## running large rasters results in the hard drive filling up with temp gribs
  ## deal with it by telling raster package to write temp files in a new defined directory within the working directory
  ## then remove the directory at the end of session
  raster_tmp_dir <- "raster_tmp"
  dir.create(raster_tmp_dir, showWarnings = F, recursive = T)
  rasterOptions(tmpdir = raster_tmp_dir)
  
  ##load imputation map of fcids
  map.fcid <- raster(file.path(results_dir, "sk_rf_2010_fcid.img", fsep = .Platform$file.sep))
  
  ##load dataframe of training metrics to attach to raster
  sk.str.attach <- read.csv(file=file.path(wkg_dir, "sk_str_attach.csv", fsep = .Platform$file.sep), head=TRUE, sep=",")
  
  ##load dataframe of observed vegetation metrics and structural attributes on training and validation plots
  ##subset to validation plots
  ##keep only select variables
  obs.training.validation <- read.csv(file=file.path(wkg_dir, "lidar_metrics_mean_training_validation.csv", fsep = .Platform$file.sep), head=TRUE,sep=",")  # averaged 3x3 plots with reference Ys built in prog1
  obs.validation <- subset(obs.training.validation, TV == "VALIDATION")
  obs.validation.2 <- obs.validation[,c("X","POLY250ID","elev_mean","elev_stddev","elev_p95","elev_cv",
                                        "percentage_first_returns_above_2m","percentage_first_returns_above_mean","mean_tree_height",
                                        "dominant_tree_height","loreys_height","basal_area","gross_stem_volume",
                                        "foliage_biomass","branch_biomass","crown_biomass","bark_biomass",
                                        "wood_biomass","stem_biomass","total_biomass")]
  ##basal area stem volumne and biomas values are kg per 25m cell, averaged across 3*3 cell plot
  ##need to convert to kg per hectare multiply by 16 (625m^2 to 10000m^2)
  obs.validation.2[,12:20] <- obs.validation.2[,12:20]*16
  ##Then divide by 1000 to convert biomass kg/ha to tons/ha for comparision to Boudewyn et al. 2007 andBeaudoin et al. 2013
  min(obs.validation.2$total_biomass/1000)
  max(obs.validation.2$total_biomass/1000)
  mean(obs.validation.2$total_biomass/1000)
  ##or divide by 1000 to convert biomass kg/ha to tons/ha for comparision to Boudewyn et al. 2007 andBeaudoin et al. 2013
  min(obs.validation.2$total_biomass/1000)
  max(obs.validation.2$total_biomass/1000)
  mean(obs.validation.2$total_biomass/1000)
  
  
  ##load shapefile of training and validation plots
  ##subset to validation models 
  ##create POLY250ID vector
  training.validation.polys <- readOGR(dsn = wkg_dir, layer = paste(zone,"_poly_250m_training_validation",sep = '')) ## 4340 polygons
  validation.polys <- subset(training.validation.polys, TV == "VALIDATION" )
  POLY250ID <- validation.polys$POLY250ID
  
  
  ##################################################
  ##extract predicted values on validation plots
  ##map of fcids is too big to attach strucutre variables to entire map on a desktop computer
  ##solution is to crop and mask fcid map one validation plot at a time, attach strucutre dataframe, extract mean values for each plot, and assign plot id
  ##do it as a for loop to get dataframe of predicted values for all validation plots

  
  ##run on one polygon to get column names for output dataframe
  map.fcid.crop <- crop(map.fcid, validation.polys[2,])
  map.fcid.clip <- mask(map.fcid.crop, map.fcid.crop)
  map.fcid.clip.sgdf <- as(map.fcid.clip, "SpatialGridDataFrame")
  names(map.fcid.clip.sgdf)[names(map.fcid.clip.sgdf)=="sk_rf_2010_fcid"] <- "FCID"
  map.ancil.stack <- suppressWarnings(stack(merge(map.fcid.clip.sgdf,sk.str.attach,by="FCID", all = FALSE)))   # merge by FCID the predicted IDs (training IDs) and the ancillary data with forest attributes   
  ##extract mean of exvars wieghted by nomralized proportion contribution of each cell in polygon
  ##wieghts within each polygon sum to 1
  extract.vars.names  <- names(extract(map.ancil.stack, validation.polys[2,], fun = mean, na.rm = FALSE, weights = FALSE, normalizeWeights = TRUE,
                                       cellnumbers = FALSE, small = FALSE, df = TRUE, factors = FALSE, sp = FALSE))
  
  ##create empty dataframe for mean variables rows equal number of polygons and columns equal to output above
  ##name columns of empty dataframe
  ##create vector of POLY250ID
  extract.vars <- data.frame(matrix(NA, nrow = length(validation.polys), ncol = length(extract.vars.names)))
  names(extract.vars) <- extract.vars.names
  POLY250ID <- validation.polys$POLY250ID
  
  ##for loop
  ##crop and mask fcid for one validation plot at a time
  ##attach structure dataframe to fcids within plot
  ##extract mean values for each plot
  # for(i in 1:length(validation.polys)){
  for (i in 2) {   # bc we only predicted on a small subregion covering validation plot nr 2
    map.fcid.crop <- crop(map.fcid, validation.polys[i,])
    map.fcid.clip <- mask(map.fcid.crop, map.fcid.crop)
    map.fcid.clip.sgdf <- as(map.fcid.clip, "SpatialGridDataFrame")
    names(map.fcid.clip.sgdf)[names(map.fcid.clip.sgdf)=="sk_rf_2010_fcid"] <- "FCID"
    map.ancil.stack <- suppressWarnings(stack(merge(map.fcid.clip.sgdf,sk.str.attach,by="FCID", all = FALSE))) # merge by FCID the predicted IDs (training IDs) and the ancillary data with forest attributes   
    ##extract mean of exvars wieghted by nomralized proportion contribution of each cell in polygon
    ##weights within each polygon sum to 1. "extract.vars" contains the validation predictions (imputations at the 25 m level upscaled at the 75 m level by averaging)
    extract.vars[i,]  <- extract(map.ancil.stack, validation.polys[i,], fun = mean, na.rm = FALSE, weights = FALSE, normalizeWeights = TRUE,
                                 cellnumbers = FALSE, small = FALSE, df = TRUE, factors = FALSE, sp = FALSE)
  } # ends i loop
  
  ##bind POLY250ID to predictions and remove uneeded ID and FCID columns
  validation.pred <- cbind(POLY250ID,extract.vars) 
  validation.pred <- validation.pred[,-c(2,3)]
  ##multiple basal area volumne and biomass values by 16 to get per hectare
  names(validation.pred)
  validation.pred[,11:19] <- validation.pred[,11:19]*16
  ##remove X column from validation observations
  ##change names of observed values except POLY250ID to end in suffix .o
  obs.validation.3 <- obs.validation.2[,-1]
  names1 <- paste(names(obs.validation.3),"o", sep = ".")
  names(obs.validation.3) <- names1
  
  ##merge obs and pred values for validation plots
  validation.obs.pred <- merge(validation.pred, obs.validation.3, by.x="POLY250ID", by.y = "POLY250ID.o")
  
  ############ TO DEL
  ## added NA subsetting bc we only tested on one single plot (NA for the rest)
  validation.obs.pred <- subset(validation.obs.pred, elev_mean != "NA") 
  
  ##export obs anf pred values for valiadtion plots
  write.csv(validation.obs.pred, file = file.path(wkg_dir, "rf_2010_validation_obs_pred.csv", fsep = .Platform$file.sep))
  
  
  #################################################################################
  ##accuracy stats function
  gmfr.stats <-
    function(x,y)
    { 
      maxrange<-max(max(x),max(y))
      minrange<-min(min(x),min(y))
      rng<-maxrange-minrange
      n<-length(x)
      ssd<-sum((x-y)^2)        				# Sum of square difference;
      msd<-ssd/n 									# Mean square difference;
      xmean<-sum(x)/n 
      ymean<-sum(y)/n								# Mean values of x and y;
      xrange<-range(x)            #range values of x and y;
      yrange<-range(y)
      xstdev<-sd(x)          #standard devaition values for x and y;
      ystdev<-sd(y)
      s_xx<-sum((x-xmean)^2)	
      s_yy<-sum((y-ymean)^2) 
      s_xy<-sum((x-xmean)*(y-ymean))
      s_xyx<-sum(((x-xmean)+(y-xmean))^2)
      d<-1-(ssd/s_xyx)								# Calculating Willmott's Index of Agreement
      rsq<-(s_xy)^2/(s_xx*s_yy)							# R-square;
      spod<-sum((abs(xmean-ymean)+abs(x-xmean))*(abs(xmean-ymean)+abs(y-ymean))) # Sum of potential difference;
      ac<-1-ssd/spod								# Agreement coefficient;
      b_yvsx<-sqrt(s_yy/s_xx)						# Estimate of b for GMFR regression y=a+bx;
      a_yvsx<-ymean-b_yvsx*xmean						# Estimate of a for GMFR regression y=a+bx;
      b_xvsy<-sqrt(s_xx/s_yy)						# Estimate of b for GMFR regression x=a+by;
      a_xvsy<-xmean-b_xvsy*ymean						# Estimate of a for GMFR regression x=a+by;
      yhat<-a_yvsx+b_yvsx*x 							# Prediction of y;  
      # lmfit_map<-lm(yhat~x)  
      lmfit_map <- NULL
      xhat<-a_xvsy+b_xvsy*y	                  			# Prediction of x;
      spd_uns<-sum(abs(x-xhat)*abs(y-yhat))     			# Unsystematic sum of product-difference;
      spd_sys<-ssd-spd_uns                   				# Systematic sum of product-difference;
      ac_uns<-1-spd_uns/spod                    			# Unsystematic agreement coefficient;
      ac_sys<-1-spd_sys/spod  		                  	# Systematic agreement coefficient;
      mpd_uns<-spd_uns/n 							# Unsystematic mean product-difference;
      mpd_sys<-spd_sys/n 							# Systematic mean product-difference;
      rmsd<-sqrt(msd)								# Root mean of square difference;
      nrmsd<-rmsd/rng               # Normalized root mean of square difference;
      rmpd_uns<-sqrt(mpd_uns)						# Unsystematic square root of mean product-difference;
      rmpd_sys<-sqrt(mpd_sys)						# Systematic square root of mean product-difference;
      mse_sys<-(sum((x-yhat)^2))/n		# Calculating Willmott's measures of agreement (systematic and unsystematic )
      mse_uns<-(sum((y-yhat)^2))/n
      mse<-mse_sys + mse_uns
      prop_uns<-mpd_uns/msd		# calculating proportion of error that's due to unsystematic differences (scatter)
      prop_sys<-mpd_sys/msd		# calculating proportion of error that's due to systematic differences (Y could be modeled from X)
      return(list(lmfit = lmfit_map, ac = list(ac = ac, uns = ac_uns,sys = ac_sys), prop = list(sys = prop_sys, uns = prop_uns),
                  rmsd = rmsd, nrmsd = nrmsd, mse = list(mse = mse, sys = mse_sys, uns = mse_uns), rsq = rsq, yhat = yhat, d = d,maxrange = maxrange,
                  xmean = xmean, ymean = ymean, xrange = xrange, yrange = yrange, xstdev = xstdev, ystdev = ystdev,
                  byvsx = b_yvsx, ayvsx = a_yvsx, bxvsy = b_xvsy, axvsy = a_xvsy))
    }
  
  #####################
  ##accuracy metrics for selected lidar vegetation metrics and derived forest strucutral attributes 
  elev_mean <- unlist(gmfr.stats(validation.obs.pred$elev_mean,validation.obs.pred$elev_mean.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  elev_stddev <- unlist(gmfr.stats(validation.obs.pred$elev_stddev,validation.obs.pred$elev_stddev.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  elev_p95 <- unlist(gmfr.stats(validation.obs.pred$elev_p95,validation.obs.pred$elev_p95.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  elev_cv <- unlist(gmfr.stats(validation.obs.pred$elev_cv,validation.obs.pred$elev_cv.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  percentage_first_returns_above_2m <- unlist(gmfr.stats(validation.obs.pred$percentage_first_returns_above_2m,validation.obs.pred$percentage_first_returns_above_2m.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  percentage_first_returns_above_mean <- unlist(gmfr.stats(validation.obs.pred$percentage_first_returns_above_mean,validation.obs.pred$percentage_first_returns_above_mean.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  mean_tree_height <- unlist(gmfr.stats(validation.obs.pred$mean_tree_height,validation.obs.pred$mean_tree_height.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  dominant_tree_height <- unlist(gmfr.stats(validation.obs.pred$dominant_tree_height,validation.obs.pred$dominant_tree_height.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  loreys_height <- unlist(gmfr.stats(validation.obs.pred$loreys_height,validation.obs.pred$loreys_height.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  basal_area <- unlist(gmfr.stats(validation.obs.pred$basal_area,validation.obs.pred$basal_area.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  gross_stem_volume <- unlist(gmfr.stats(validation.obs.pred$gross_stem_volume,validation.obs.pred$gross_stem_volume.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  foliage_biomass  <- unlist(gmfr.stats(validation.obs.pred$foliage_biomass,validation.obs.pred$foliage_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  branch_biomass  <- unlist(gmfr.stats(validation.obs.pred$branch_biomass,validation.obs.pred$branch_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  crown_biomass <- unlist(gmfr.stats(validation.obs.pred$crown_biomass,validation.obs.pred$crown_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  bark_biomass <- unlist(gmfr.stats(validation.obs.pred$bark_biomass,validation.obs.pred$bark_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  wood_biomass <- unlist(gmfr.stats(validation.obs.pred$wood_biomass,validation.obs.pred$wood_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  stem_biomass <- unlist(gmfr.stats(validation.obs.pred$stem_biomass,validation.obs.pred$stem_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  total_biomass <- unlist(gmfr.stats(validation.obs.pred$total_biomass,validation.obs.pred$total_biomass.o)[c("ac","nrmsd","rsq","rmsd","xmean","ymean","xrange","yrange","xstdev","ystdev","byvsx","ayvsx","bxvsy","axvsy")])
  
  
  ## merge accuract metrics for all lidar vegetation metrics and deirved strucutral atttru=ibutes
  ## merge the ac rmsd rsq and ks statistics from the gmfr and ecdf functions
  ## make seperate living dead objects
  gmfr.stats.str <- rbind(elev_mean,elev_stddev,elev_p95,elev_cv,percentage_first_returns_above_2m,percentage_first_returns_above_mean,
                          mean_tree_height,dominant_tree_height,loreys_height,basal_area,gross_stem_volume,
                          foliage_biomass,branch_biomass,crown_biomass,bark_biomass,wood_biomass,stem_biomass,total_biomass)
  
  
  write.csv(gmfr.stats.str, file = file.path(wkg_dir, "rf_2010_aa_metrics.csv", fsep = .Platform$file.sep))
  
  ##delete temp folder with temp rasters
  unlink(raster_tmp_dir, recursive = T, force = T)
  
  
  ##############################################
  ##plot observed versus predicted values for 6 lidar metrics
  ##lots of plots makes it desirable for semi-tranparency to using ggplot
  ##scatter plots of obs versus pred
  ##with a 1:1 line in black
  ##with gmfr line in red
  ##call multiplot function
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
#   ##read obs and pred values for valiadtion plots
#   validation.obs.pred <- read.csv(file = "J:/sk_models/v1/accass_2010/sk_rf_2010_validation_obs_pred.csv",head=TRUE,sep=",")
#   ##read gmfr.stats.str
#   gmfr.stats.str <- read.csv(file = "J:/sk_models/v1/accass_2010/sk_rf_2010_aa_metrics.csv",head=TRUE,sep=",")
#   
  ##ggplot of 6 lidar vars
  p1 <- ggplot(validation.obs.pred, aes(elev_mean.o, elev_mean))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="mean height (m)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[1,19], slope = gmfr.stats.str[1,18], colour = 'red')+
    scale_x_continuous(limits=c(2,17))+
    scale_y_continuous(limits=c(2,17))+
    geom_text(data = NULL, x = 9.5, y = 17, label = "r^2 == 0.52", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA))
  
  p2 <- ggplot(validation.obs.pred, aes(elev_stddev.o, elev_stddev))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="stdev of height (m)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[2,19], slope = gmfr.stats.str[2,18], colour = 'red')+
    scale_x_continuous(limits=c(0,8))+
    scale_y_continuous(limits=c(0,8))+
    geom_text(data = NULL, x = 4, y = 8, label = "r^2 == 0.48", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA))
  
  p3 <- ggplot(validation.obs.pred, aes(elev_p95.o, elev_p95))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="95th percentile of height (m)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[3,19], slope = gmfr.stats.str[3,18], colour = 'red')+
    scale_x_continuous(limits=c(3,26))+
    scale_y_continuous(limits=c(3,26))+
    geom_text(data = NULL, x = 14, y = 26, label = "r^2 == 0.54", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA))
  
  p4 <- ggplot(validation.obs.pred, aes(elev_cv.o, elev_cv))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="CV of height (%)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[4,19], slope = gmfr.stats.str[4,18], colour = 'red')+
    scale_x_continuous(limits=c(0.1,0.85))+
    scale_y_continuous(limits=c(0.1,0.85))+
    geom_text(data = NULL, x = 0.475, y = 0.85, label = "r^2 == 0.42", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA))
  
  p5 <- ggplot(validation.obs.pred, aes(percentage_first_returns_above_2m.o, percentage_first_returns_above_2m))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="cover above 2m (%)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[5,19], slope = gmfr.stats.str[5,18], colour = 'red')+
    scale_x_continuous(limits=c(0,100))+
    scale_y_continuous(limits=c(0,100))+
    geom_text(data = NULL, x = 50, y = 100, label = "r^2 == 0.69", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p6 <- ggplot(validation.obs.pred, aes(percentage_first_returns_above_mean.o, percentage_first_returns_above_mean))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="cover above mean height (%)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[6,19], slope = gmfr.stats.str[6,18], colour = 'red')+
    scale_x_continuous(limits=c(0,60))+
    scale_y_continuous(limits=c(0,60))+
    geom_text(data = NULL, x = 30, y = 60, label = "r^2 == 0.68", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  ## ggplot 12 derived vars
  p7 <- ggplot(validation.obs.pred, aes(mean_tree_height.o, mean_tree_height))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="mean height (m)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[7,19], slope = gmfr.stats.str[7,18], colour = 'red')+
    scale_x_continuous(limits=c(6,18))+
    scale_y_continuous(limits=c(6,18))+
    geom_text(data = NULL, x = 12, y = 18, label = "r^2 == 0.56", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p8 <- ggplot(validation.obs.pred, aes(dominant_tree_height.o, dominant_tree_height))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="dominant height (m)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[8,19], slope = gmfr.stats.str[8,18], colour = 'red')+
    scale_x_continuous(limits=c(6,28))+
    scale_y_continuous(limits=c(6,28))+
    geom_text(data = NULL, x = 17, y = 28, label = "r^2 == 0.58", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p9 <- ggplot(validation.obs.pred, aes(loreys_height.o, loreys_height))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title="lorey's height (m)")+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[9,19], slope = gmfr.stats.str[9,18], colour = 'red')+
    scale_x_continuous(limits=c(4,22))+
    scale_y_continuous(limits=c(4,22))+
    geom_text(data = NULL, x = 12.5, y = 22, label = "r^2 == 0.55", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p10 <- ggplot(validation.obs.pred, aes(basal_area.o, basal_area))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("basal"," ","area"," ","(",m^2,"/"," ","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[10,19], slope = gmfr.stats.str[10,18], colour = 'red')+
    scale_x_continuous(limits=c(0,50))+
    scale_y_continuous(limits=c(0,50))+
    geom_text(data = NULL, x = 25, y = 50, label = "r^2 == 0.57", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p11 <- ggplot(validation.obs.pred, aes(gross_stem_volume.o, gross_stem_volume))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("gross"," ","stem"," ","volume"," ","(",m^3,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[11,19], slope = gmfr.stats.str[11,18], colour = 'red')+
    scale_x_continuous(limits=c(0,250))+
    scale_y_continuous(limits=c(0,250))+
    geom_text(data = NULL, x = 125, y = 250, label = "r^2 == 0.53", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p12 <- ggplot(validation.obs.pred, aes(foliage_biomass.o, foliage_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("foliage"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[12,19], slope = gmfr.stats.str[12,18], colour = 'red')+
    scale_x_continuous(limits=c(0,7000))+
    scale_y_continuous(limits=c(0,7000))+
    geom_text(data = NULL, x = 3500, y = 7000, label = "r^2 == 0.65", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p13 <- ggplot(validation.obs.pred, aes(branch_biomass.o, branch_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("branch"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[13,19], slope = gmfr.stats.str[13,18], colour = 'red')+
    scale_x_continuous(limits=c(0,32000))+
    scale_y_continuous(limits=c(0,32000))+
    geom_text(data = NULL, x = 16000, y = 32000, label = "r^2 == 0.57", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  
  p14 <- ggplot(validation.obs.pred, aes(crown_biomass.o, crown_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("crown"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[14,19], slope = gmfr.stats.str[14,18], colour = 'red')+
    scale_x_continuous(limits=c(0,40000))+
    scale_y_continuous(limits=c(0,40000))+
    geom_text(data = NULL, x = 20000, y = 40000, label = "r^2 == 0.59", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p15 <- ggplot(validation.obs.pred, aes(bark_biomass.o, bark_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("bark"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[15,19], slope = gmfr.stats.str[15,18], colour = 'red')+
    scale_x_continuous(limits=c(0,24000))+
    scale_y_continuous(limits=c(0,24000))+
    geom_text(data = NULL, x = 12000, y = 24000, label = "r^2 == 0.54", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p16 <- ggplot(validation.obs.pred, aes(wood_biomass.o, wood_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("wood"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[16,19], slope = gmfr.stats.str[16,18], colour = 'red')+
    scale_x_continuous(limits=c(0,192000))+
    scale_y_continuous(limits=c(0,192000))+
    geom_text(data = NULL, x = 96000, y = 192000, label = "r^2 == 0.52", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p17 <- ggplot(validation.obs.pred, aes(stem_biomass.o, stem_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("stem"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[17,19], slope = gmfr.stats.str[17,18], colour = 'red')+
    scale_x_continuous(limits=c(0,208000))+
    scale_y_continuous(limits=c(0,208000))+
    geom_text(data = NULL, x = 104000, y = 208000, label = "r^2 == 0.52", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  p18 <- ggplot(validation.obs.pred, aes(total_biomass.o, total_biomass))+
    geom_point(colour = 'gray10', size = 3, alpha = 0.2)+
    labs(x="observed", y="predicted", title=expression(paste("total"," ","biomass"," ","(",kg,"/","ha",")")))+
    geom_abline(intercept = 0, slope = 1, colour = 'black')+
    geom_abline(intercept = gmfr.stats.str[18,19], slope = gmfr.stats.str[18,18], colour = 'red')+
    scale_x_continuous(limits=c(0,560000))+
    scale_y_continuous(limits=c(0,560000))+
    geom_text(data = NULL, x = 280000, y = 560000, label = "r^2 == 0.50", parse = TRUE)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA))
  
  
  ## multiple plot graph of 6 lidar vars plus 3 ancillary derived attributes
  ## multiplot order is wierd up to down left to right
  multiplot(p1, p3, p10, p2, p5, p11, p4, p6, p18, cols=3)
  ## export as width 800 ht 800
  
  
  ## multiple plot graph of 9 remaining ancillary derived attributes
  ## multiplot order is wierd up to down left to right
  multiplot(p7, p12, p15, p8, p13, p16, p9, p14, p17, cols=3)
  ## export as width 800 ht 800
  
  
  
  
  
  plot(validation.obs.pred$total_biomass.o,validation.obs.pred$total_biomass, pch = 16,
       col = "gray", main = "total biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,560000), ylim=c(0,560000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[18,18],b=gmfr.stats.str[18,17], col = "red")
  text(280000,560000, "R2 = 0.50")
  
  
  
  
  
  
  
  ##############################################
  ##plot observed versus predicted for 6 lidar metrics
  ##within a 1:1 line in black
  ##with gmfr line in red
  
  par(mfrow=c(2,3))
  plot(validation.obs.pred$elev_mean.o,validation.obs.pred$elev_mean, pch = 16, col = "gray", main = "mean height (m)",
       xlab = "observed", ylab = "predicted", xlim = c(2,17), ylim=c(2,17))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[1,19],b=gmfr.stats.str[1,18], col = "red")
  text(10,17, "R2 = 0.52")
  
  plot(validation.obs.pred$elev_stddev.o,validation.obs.pred$elev_stddev, pch = 16, col = "gray", main = "stdev of height (m)",
       xlab = "observed", ylab = "predicted", xlim = c(0,8), ylim=c(0,8))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[2,18],b=gmfr.stats.str[2,17], col = "red")
  text(4,8, "R2 = 0.48")
  
  plot(validation.obs.pred$elev_p95.o,validation.obs.pred$elev_p95, pch = 16, col = "gray", main = "95th percentile height (m)",
       xlab = "observed", ylab = "predicted", xlim = c(0,26), ylim=c(0,26))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[3,18],b=gmfr.stats.str[3,17], col = "red")
  text(12.5,26, "R2 = 0.54")
  
  plot(validation.obs.pred$elev_cv.o,validation.obs.pred$elev_cv, pch = 16, col = "gray", main = "CV of height (percent)",
       xlab = "observed", ylab = "predicted", xlim = c(0,0.85), ylim=c(0,0.85))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[4,18],b=gmfr.stats.str[4,17], col = "red")
  text(0.4,0.85, "R2 = 0.42")
  
  plot(validation.obs.pred$percentage_first_returns_above_2m.o,validation.obs.pred$percentage_first_returns_above_2m, 
       pch = 16, col = "gray", main = "cover above 2m (%)", xlab = "observed", ylab = "predicted", xlim = c(0,100), ylim=c(0,100))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[5,18],b=gmfr.stats.str[5,17], col = "red")
  text(50,100, "R2 = 0.69")
  
  plot(validation.obs.pred$percentage_first_returns_above_mean.o,validation.obs.pred$percentage_first_returns_above_mean, 
       pch = 16, col = "gray", main = "cover above mean height (%)", xlab = "observed", ylab = "predicted", xlim = c(0,60), ylim=c(0,60))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[6,18],b=gmfr.stats.str[6,17], col = "red")
  text(30,60, "R2 = 0.68")
  
  ##############################################
  ##plot observed versus predicted for 12 derived forest strucutre 6 lidar metrics
  ##within a 1:1 line in black
  ##with gmfr line in red
  par(mfrow=c(4,3))
  plot(validation.obs.pred$mean_tree_height.o,validation.obs.pred$mean_tree_height, pch = 16,
       col = "gray", main = "mean height (m)", xlab = "observed", ylab = "predicted", xlim = c(6,17), ylim=c(6,17))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[7,18],b=gmfr.stats.str[7,17], col = "red")
  text(11,17, "R2 = 0.56")
  
  plot(validation.obs.pred$dominant_tree_height.o,validation.obs.pred$dominant_tree_height, pch = 16,
       col = "gray", main = "dominant height (m)", xlab = "observed", ylab = "predicted", xlim = c(5,28), ylim=c(5,28))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[8,18],b=gmfr.stats.str[8,17], col = "red")
  text(15,28, "R2 = 0.58")
  
  plot(validation.obs.pred$loreys_height.o,validation.obs.pred$loreys_height, pch = 16,
       col = "gray", main = "lorey's height (m)", xlab = "observed", ylab = "predicted", xlim = c(4,22), ylim=c(4,22))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[9,18],b=gmfr.stats.str[9,17], col = "red")
  text(12.5,22, "R2 = 0.55")
  
  plot(validation.obs.pred$basal_area.o,validation.obs.pred$basal_area, pch = 16,
       col = "gray", main = "basal area (m^2/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,50), ylim=c(0,50))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[10,18],b=gmfr.stats.str[10,17], col = "red")
  text(25,50, "R2 = 0.57")
  
  plot(validation.obs.pred$gross_stem_volume.o,validation.obs.pred$gross_stem_volume, pch = 16,
       col = "gray", main = "gross stem volume (m^3/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,250), ylim=c(0,250))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[11,18],b=gmfr.stats.str[11,17], col = "red")
  text(125,250, "R2 = 0.53")
  
  plot(validation.obs.pred$foliage_biomass.o,validation.obs.pred$foliage_biomass, pch = 16,
       col = "gray", main = "foliage biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,7000), ylim=c(0,7000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[12,18],b=gmfr.stats.str[12,17], col = "red")
  text(3500,7000, "R2 = 0.65")
  
  plot(validation.obs.pred$branch_biomass.o,validation.obs.pred$branch_biomass, pch = 16,
       col = "gray", main = "branch biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,32000), ylim=c(0,32000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[13,18],b=gmfr.stats.str[13,17], col = "red")
  text(16000,32000, "R2 = 0.57")
  
  plot(validation.obs.pred$crown_biomass.o,validation.obs.pred$crown_biomass, pch = 16,
       col = "gray", main = "crown biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,40000), ylim=c(0,40000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[14,18],b=gmfr.stats.str[14,17], col = "red")
  text(20000,40000, "R2 = 0.59")
  
  plot(validation.obs.pred$bark_biomass.o,validation.obs.pred$bark_biomass, pch = 16,
       col = "gray", main = "bark biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,24000), ylim=c(0,24000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[15,18],b=gmfr.stats.str[15,17], col = "red")
  text(12000,24000, "R2 = 0.54")
  
  plot(validation.obs.pred$wood_biomass.o,validation.obs.pred$wood_biomass, pch = 16,
       col = "gray", main = "wood biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,192000), ylim=c(0,192000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[16,18],b=gmfr.stats.str[16,17], col = "red")
  text(96000,192000, "R2 = 0.52")
  
  plot(validation.obs.pred$stem_biomass.o,validation.obs.pred$stem_biomass, pch = 16,
       col = "gray", main = "stem biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,208000), ylim=c(0,208000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[17,18],b=gmfr.stats.str[17,17], col = "red")
  text(104000,208000, "R2 = 0.52")
  
  plot(validation.obs.pred$total_biomass.o,validation.obs.pred$total_biomass, pch = 16,
       col = "gray", main = "total biomass (kg/ha)", xlab = "observed", ylab = "predicted", xlim = c(0,560000), ylim=c(0,560000))
  abline(a=0,b=1, col = "black")
  abline(a=gmfr.stats.str[18,18],b=gmfr.stats.str[18,17], col = "red")
  text(280000,560000, "R2 = 0.50")
  
  # clock UTM zone time
  temp.toc <- proc.time()-temp.tic[3]
  print(paste(zone,"elapsed time:",seconds_to_period(temp.toc[3])))
  
}

# stopCluster(cl)

# clock global time
toc <- proc.time()-tic[3]
print(paste("Total elapsed time:",seconds_to_period(toc[3])))


