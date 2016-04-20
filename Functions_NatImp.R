##########################################################
# Project Name: CA_IMPUTATION
# Author: Giona Matasci          
# File Name: Functions_NatImp.R                             
##########################################################

## function to copy/paste UTM zone data (NTEMS folder structure) from external HDD to a local directory
copy_paste_UTMdata <- function(disk.names, zone, cng.var.names.long, local.dir){
   
  tic.copy.UTMdata <- proc.time()
  
  list.missing.files <- vector("list", 0)
  zone.nr <- substr(zone, 4, nchar(zone))
  
  for (dn in disk.names) {
   list.dir.in.disk <- list.dirs(path = dn, full.names = FALSE, recursive = FALSE)
   
   ## checks if our zone is in the disk, if not true goes to next round of the loop with "next" command
   if ( !zone %in% gsub('_', '', list.dir.in.disk) ) {
     next  
   }
  
   #### copy raw BAP Landsat file for 2010 to compute TC
   cmd <- sprintf('relative.path <- file.path("UTM_%s", "Results", "proxy_values", fsep = .Platform$file.sep)', zone.nr)  
   eval(parse(text=cmd))
   target.dir <- file.path(local.dir, relative.path, fsep = .Platform$file.sep)
   dir.create(target.dir, showWarnings = F, recursive = T)
   
   cmd <- sprintf('file.to.copy <- file.path("%s", relative.path, "SRef_UTM%s_2010_proxy_v2", fsep = .Platform$file.sep)', dn, zone.nr)
   eval(parse(text=cmd))
   file.to.copy.hdr <- paste (file.to.copy, '.hdr', sep = "")
   file.to.copy.dat <- paste (file.to.copy, '.dat', sep = "")
   files.to.copy <- list(file.to.copy.hdr, file.to.copy.dat)
   ##----- to test script rapidly ----------
   # files.to.copy <- list(file.to.copy.hdr)
   ##---------------------------------------
   
   if ( all(file.exists(unlist(files.to.copy))) ) {
     file.copy(from=files.to.copy, to=target.dir, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
   } else {
     warning(paste ("Missing file:", file.to.copy, sep = " "))
     list.missing.files <- append(list.missing.files, file.to.copy)
   }
   
   #### copy change metrics files
   cmd <- sprintf('relative.path <- file.path("UTM_%s", "Results", "Change_metrics", fsep = .Platform$file.sep)', zone.nr)  
   eval(parse(text=cmd))
   target.dir <- file.path(local.dir, relative.path, fsep = .Platform$file.sep)
   dir.create(target.dir, showWarnings = F, recursive = T)
  
   for (var in cng.var.names.long) {
     cmd <- sprintf('file.to.copy <- file.path("%s", relative.path, "SRef_%s_%s", fsep = .Platform$file.sep)', dn, zone.nr, var)
     eval(parse(text=cmd))
     
     file.to.copy.hdr <- paste (file.to.copy, '.hdr', sep = "")
     file.to.copy.dat <- paste (file.to.copy, '.dat', sep = "")
     files.to.copy <- list(file.to.copy.hdr, file.to.copy.dat)
     ##----- to test script rapidly ----------
     # files.to.copy <- list(file.to.copy.hdr)
     ##---------------------------------------
     
     if ( all(file.exists(unlist(files.to.copy))) ) {
       file.copy(from=files.to.copy, to=target.dir, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
     } else {
       warning(paste ("Missing file:", file.to.copy, sep = " "))
       list.missing.files <- append(list.missing.files, file.to.copy)
     }
     
   }
   
   #### copy change attribution file
   cmd <- sprintf('relative.path <- file.path("UTM_%s", "Results", "change_attribution", fsep = .Platform$file.sep)', zone.nr)  
   eval(parse(text=cmd))
   target.dir <- file.path(local.dir, relative.path, fsep = .Platform$file.sep)
   dir.create(target.dir, showWarnings = F, recursive = T)
   
   cmd <- sprintf('file.to.copy <- file.path("%s", relative.path, "UTM%s_Change_attribution_complete", fsep = .Platform$file.sep)', dn, zone.nr)
   eval(parse(text=cmd))
   file.to.copy.hdr <- paste (file.to.copy, '.hdr', sep = "")
   file.to.copy.dat <- paste (file.to.copy, '.dat', sep = "")
   files.to.copy <- list(file.to.copy.hdr, file.to.copy.dat)
   ##----- to test script rapidly ----------
   # files.to.copy <- list(file.to.copy.hdr)
   ##---------------------------------------
   
   if ( all(file.exists(unlist(files.to.copy))) ) {
     file.copy(from=files.to.copy, to=target.dir, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
   } else {
     warning(paste ("Missing file:", file.to.copy, sep = " "))
     list.missing.files <- append(list.missing.files, file.to.copy)
   }
  
  }
  
  toc <- proc.time()-tic.copy.UTMdata[3]
     
  return( list(list.missing.files=list.missing.files, toc=toc) )
   
}
 
## function to convert Landsat band values to Tasseled Cap components values
bands_2_TC <- function(bands, sens='ETM+'){
   
   if(sens=='ETM+') {
     ## TC transform (with ETM+ coeff)
     tc.coef <- t(matrix(c(
       # Tasseled cap coefficients for Landsat 7 ETM+ at-satellite reflectance from HWY+2002
       # Band 1     Band 2       Band 3     Band 4     Band 5        Band 7     Index
       0.3561,     0.3972,      0.3904,    0.6966,    0.2286,       0.1596,    #  Brightness       
       -0.3344,    -0.3544,     -0.4556,    0.6966,   -0.0242,      -0.2630,    #  Greenness       
       0.2626,     0.2141,      0.0926,    0.0656,   -0.7629,      -0.5388    #  Wetness          
     ), ncol=6, byrow=TRUE))
   } else if (sens=='TM') {
     ## TC transform (with TM coeff)
     tc.coef <- t(matrix(c(
       # TM Tasseled Cap Equivalent Transformation Matrix for Band Reflectance Factor from Crist1985
       # Band 1     Band 2       Band 3     Band 4     Band 5        Band 7     Index
       0.2043,     0.4158,      0.5524,    0.5741,    0.3124,       0.2303,    #  Brightness
       -0.1603,    -0.2819,     -0.4934,    0.7940,    0.0002,      -0.1446,    #  Greenness
       0.0315,     0.2021,      0.3102,    0.1594,    0.6806,      -0.6109    #  Wetness          
     ), ncol=6, byrow=TRUE))
   }
   
   TCs <- bands %*% tc.coef
   TCang <- atan2(TCs[,2], TCs[,1])
   TCdist <- sqrt(TCs[,2]**2 + TCs[,1]**2)
   
   return( list(TCs=TCs, TCang=TCang, TCdist=TCdist) )
   
}

## function to compute the hours of day light in a certain day of the year as a function of latitude
lat_2_daylight <- function(lat, j){

  ## 3 examples
#   jul.day.vect <- 1:365
#   lat <- 40
#   daylight <- lat_2_daylight(lat, jul.day.vect) 
#   plot(jul.day.vect, daylight)
#   
#   jul.day <- 172
#   lat.vect <- 0:70
#   daylight <- lat_2_daylight(lat.vect, jul.day) 
#   plot(lat.vect, daylight)
  
#   jul.day.vect <- c(21,	52,	80,	111,	141,	172,	202,	233,	264,	294,	325,	355)
#   lat.vect <- 0:70
#   daylight.tab <- array(NaN, c(length(lat.vect), length(jul.day.vect)))
#   jidx <- 1
#   for (jul.day in jul.day.vect) {
#     daylight.tab[, jidx] <- lat_2_daylight(lat.vect, jul.day)
#     jidx <- jidx+1
#   }
#   mean.daylight <- rowMeans(daylight.tab)
#   plot(mean.daylight)
  
  P = asin( .39795*cos( .2163108 + 2*atan( .9671396*tan(.00860*(j-186)) ) ) )
  D = 24 - (24/pi)*acos( ( sin(0.8333*pi/180) + sin(lat*pi/180)*sin(P) ) / ( cos(lat*pi/180)*cos(P) ) )
  return (D)
}

## function to investigate the behavior of Random Forest (regression) parameters
perf_oob_Random_Forest <- function(X, Y, params, seed){
  if (hasArg(seed)) {   # check if argument has been specified
    set.seed(seed)
  }
  rmsd <- array(rep(NaN, length(params$nrtrees.vect)*length(params$mtry.vect)*length(params$nodesize.vect)), c(length(params$nrtrees.vect), length(params$mtry.vect), length(params$nodesize.vect)))
  rsq <- rmsd
  i <- 1
  for (nrtrees in params$nrtrees.vect) {
    print(sprintf('Nr. trees: %i [%i/%i], max: %i', nrtrees, i, length(params$nrtrees.vect), params$nrtrees.vect[length(params$nrtrees.vect)]))
    j <-1
    for (mtries in params$mtry.vect) {
      l <- 1
      for (nodesizes in params$nodesize.vect) {
        rf.model <- randomForest(x=X, y=Y, ntree=nrtrees, mtry=mtries, nodesize=nodesizes)
        rmsd[i,j,l] <- sqrt(rf.model$mse[length(rf.model$mse)])  # last element of rf.model$mse is the mse (or rsq) of the whole forest 
        rsq[i,j,l] <- rf.model$rsq[length(rf.model$rsq)] 
        l <- l+1
      }
      j <- j+1
    }
    i <- i+1
  }
  return(list(rmsd=rmsd, rsq=rsq))
}

## function to assess the performance of a regression model
regr_metrics <- function(x,y) { 
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
    lmfit_map<-lm(yhat~x)  
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
    return(
      list(
        lmfit = lmfit_map, 
        ac = list(ac = ac, uns = ac_uns,sys = ac_sys), 
        prop = list(sys = prop_sys, uns = prop_uns),
        rmsd = rmsd, 
        nrmsd = nrmsd, 
        mse = list(mse = mse, sys = mse_sys, uns = mse_uns), 
        rsq = rsq, 
        yhat = yhat, 
        d = d,
        maxrange = maxrange,
        xmean = xmean, 
        ymean = ymean, 
        xrange = xrange,
        yrange = yrange, 
        xstdev = xstdev, 
        ystdev = ystdev,
        byvsx = b_yvsx, 
        ayvsx = a_yvsx, 
        bxvsy = b_xvsy, 
        axvsy = a_xvsy
      )
    )
}

## function to plot a matrix as an image with axis tick labels and colorbar (from http://www.phaget4.org/R/image_matrix.html) 
my_image_plot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}

## function to convert row names of dataframes to FCID
RownamesToFCID <- function(indataframe){
  rownames(indataframe) <- indataframe[,"FCID"]
  indataframe <- indataframe[,2:ncol(indataframe)]
  return(indataframe)
}

