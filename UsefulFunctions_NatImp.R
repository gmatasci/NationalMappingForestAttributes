##########################################################
# Project Name: CA_IMPUTATION
# Author: Giona Matasci          
# File Name: UsefulFunctions_NatImp.R                             
# Objective: Useful functions
##########################################################

copy_paste_UTMdata <- function(disk.names, zone, cng.var.names.long, local.dir){
   
   tic.copy.UTMdata <- proc.time()
   
   #    zone.numbers <- c(7:22)
   #    zone.latitudes <- c('S', 'N', 'A') 
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
     
     cmd <- sprintf('file.to.copy <- file.path("%s", relative.path, "SRef_UTM%s_2010_proxy", fsep = .Platform$file.sep)', dn, zone.nr)
     eval(parse(text=cmd))
     file.to.copy.hdr <- paste (file.to.copy, '.hdr', sep = "")
     file.to.copy.dat <- paste (file.to.copy, '.dat', sep = "")
     files.to.copy <- list(file.to.copy.hdr, file.to.copy.dat)
     # files.to.copy <- list(file.to.copy.hdr)
     
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
       # files.to.copy <- list(file.to.copy.hdr)
       
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
     # files.to.copy <- list(file.to.copy.hdr)
     
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


