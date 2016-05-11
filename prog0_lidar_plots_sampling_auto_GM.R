##----------------------
## TO DO
##----------------------

#### STILL TO DO ####
# -DOUG: check validity of Ecozone field

#### SOLVED ####
# -V check IDs who relates to who -- now they all relate to "lidar.2.fc.dist.undist"
# -V check shapefiles produced -- FCIDs all match and relate to "lidar.2.fc.dist.undist"
# -V duplicate plot IDs (eg, in "9S_points.shp": T20_095D01_333871 and T24_095D01_333871) -- now filtered in the new version of files Geordie sent ("LOP_attr.zip" folder)
# -V duplicate plot coords -- see impact in script and just filter them myself (after "lidar <-  readOGR()" )
# -V check projections -- inserted if-else-stop to check
# -V check usage of new shp in LOP_attributes as they have more columns/points and new Ecozone column -- no errors are given
# -V check code step by step with new file names and new split of layers poly_sampling/filtering
# -V problem with FCID <- polys.250m.9pts.df@data$id when building sampling polygons: id is not matching with FCID later on (prog1) -- now it is: FCID <- lidar.2.fc.dist.undist.polys.250m.centerpts@data$FCID, so it refers to the same centerpoints file


##----------------------
## READS  
##----------------------
# - "<UTMzone>_points.shp": point shp with LiDAR 25x25 m plots
# - "<UTMzone>_trsct.shp": polyline shp with center of transect line

##----------------------
## WRITES
##----------------------
# - "<UTMzone>_lop_wkg_fc.shp": point shp with plots surrounded by forest in a 3x3 neighborhood
# - "<UTMzone>_lop_wkg_fc_dist_undist.shp": point shp with plots surrounded by forest and all either disturbed or undisturbed in a 3x3 neighborhood (MAIN OUTPUT WHOSE FCIDs REFER TO) 
# - "<UTMzone>_pt_centerpt.shp": (for prog1) point shp with center of hexagons used for sampling, after subsetting wrt forest/dist/undist (MAIN OUTPUT USED TO PRODUCE SHP "<UTMzone>_cpt_poly_250m_training_validation" IN THE NEXT STEP)
# - "<UTMzone>_poly_sampling.shp": (for prog1) polygon shp with sampling polygons (3x3 or 1x1) covering plots' squares, after subsetting wrt forest/dist/undist
# - "<UTMzone>_pt_9plots_filtering.shp": (for prog1) point shp with the LiDAR plots (3x3 = 9 plots) used for filtering

#-----------------------------------------------------------------
#-------------------------     START     -------------------------
#-----------------------------------------------------------------

print('Prog0: defining sampling grid')

rm(list=ls()) ## clear all variables

##------------------------
## LOAD GLOBAL PARAMETERS
##------------------------
param_file = "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

##!!!!!!!!!!!!!!!!!!!!
# paramsGL$zones <- c('UTM12N')
##!!!!!!!!!!!!!!!!!!!!

##----------------------------
## SCRIPT SPECIFIC PARAMETERS
##----------------------------
params0 <- list()
params0$dist.center.trans <- 300   ## threshold on distance from center of transect to eliminate plots whose acquisition was too skewed
# params0$dist.center.trans <- 100   ## smaller test threshold to speed up tests
params0$hexag.cell.size <- 250   ## spacing of the hexagonal sampling cells
# params0$hexag.cell.size <- 3500  ## larger test spacing to speed up tests

params0$min.snap.dist <- 25  ## limit distance used in function to allocate coordinates of sample points those of the nearest lidar plots

param_file_prog0 = file.path(base_wkg_dir, 'AllUTMzones_params0.Rdata', fsep = .Platform$file.sep) 
save(params0, file = param_file_prog0)


nr.clusters = length(paramsGL$zones)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end


##----------------------------
## LOAD PACKAGES
##----------------------------
list.of.packages <- c("rgdal",
                      "sp",
                      "spdep",
                      "spatstat",
                      "rgeos",
                      "maptools", 
                      "lubridate", 
                      "doParallel", 
                      "foreach"
                      )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   # named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)  ## to allow calling library() dynamically
}


nr.samples.per.zone <- array(NaN, length(paramsGL$zones))

tic <- proc.time() ## start clocking global time

cl <- makeCluster(nr.clusters)
registerDoParallel(cl)
## when using foreach the last unassigned variable can be stored in a growing dataframe (the "nr.samples.per.zone" to which the foreach result is assigned here)
nr.samples.per.zone <- foreach (z = 1:length(paramsGL$zones), .combine='rbind', .packages=list.of.packages) %dopar% {   #add .verbose=TRUE for more info when debugging
  ##---- WHEN NOT USING FOREACH, TESTING PHASE ----
  # for (z in 1:length(paramsGL$zones)) {   ## to be used when debugging (using foreach breakpoints don't work and line numbers of errors is not shown)
  ##---- WHEN NOT USING FOREACH, TESTING PHASE ----
  
  zone = paramsGL$zones[z]
  wkg_dir = file.path(base_wkg_dir,zone, fsep = .Platform$file.sep)

  if (! file.exists(wkg_dir)){dir.create(wkg_dir, showWarnings = F, recursive = T)}

  ## load lidar plot points shapefiles
  pnt_dir <- file.path(LOP_dir,'LOP_attributed', fsep = .Platform$file.sep) 
  lidar.raw <-  readOGR(dsn = pnt_dir, layer = paste(zone,"_points", sep='')) ## 1667288 lidar plots, "lidar" is of type SpatialPointsDataFrame that is an object of the S4 object system, readOGR belongs to rgdal: dsn = data source name (folder), layer = filename without extension
  
  ## filter duplicate coordinates
  duplicIndic <- duplicated(lidar.raw@coords, incomparables = FALSE)
  lidar <- lidar.raw[!duplicIndic, ]  ## lidar now contains only one point at each coordinate pair
  rm(lidar.raw)
  rm(duplicIndic)
  
  ## also load the lidar transect polyline shapefile
  trsct_dir <- file.path(LOP_dir,'LOP_transects', fsep = .Platform$file.sep) 
  lidar.transect <- readOGR(dsn = trsct_dir, layer = paste(zone, "_trsct",sep=''))
  
  ## check if coordinate systems are the same between lidar plots point shapefile and transect line shapefile
  if ( !identical(proj4string(lidar), proj4string(lidar.transect)) ) {
    stop(sprintf("%s: projections of LiDAR points shapefile and LiDAR transect shapefile loaded in prog0 do not match.", zone))
  }
  
  ## sample based on distance to edge of lidar transect
  ## keep lidar points that are no farther than params0$dist.center.trans meters (set to 300m) from the center of lidar transect
  lidar.2 <- subset(lidar, NEAR_DIST <= params0$dist.center.trans)  ## NEAR_DIST is an attribute of the layer indicating the distance from the center of the transect

  ## subset to keep only points in FC DISTURB and UNDISTURB condition classes
  lidar.2.fc <- subset(lidar.2, for_sum == 9) 
  writeOGR(lidar.2.fc , wkg_dir, paste(zone,"_lop_wkg_fc", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## save it for future checks
  lidar.2.fc.dist.undist <- subset(lidar.2.fc, dist_sum == 9 | dist_sum == 0 )
  lidar.2.fc.dist.undist@data$FCID  <- seq(1, length(lidar.2.fc.dist.undist), by=1)  # add unique identifier FCID for each of the points selected
  writeOGR(lidar.2.fc.dist.undist, wkg_dir, paste(zone,"_lop_wkg_fc_dist_undist", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)   # shapefile with points respecting all the conditions (totally disturbed or totally undisturbed forest regions)
   
  ## subset coordinates by defined minimum distance between plots
  ## sample points along hexagon lattice with 250m spacing
  ## restrict sampling to within 400 m buffer of either side of lidar transect center 
  transect.buffer <- gBuffer(lidar.transect, width=400)
  set.seed(2010)  ## spssample has a random component in where it places the grid, so needs to be initialized the same each time
  HexPts.250m <- spsample(transect.buffer, type = "hexagonal", cellsize = params0$hexag.cell.size)    ## cell size is set to 250m, that is Harold's final value

  ## generate hexagon unique IDs
  Hex250ID <- as.data.frame(seq(1, length(HexPts.250m), by=1))
  names(Hex250ID) <- "Hex250ID"  ## assign name to the column of the data frame
  HexPts.250m.SPDF <- SpatialPointsDataFrame(coordinates(HexPts.250m), data=Hex250ID, proj4string=CRS(proj4string(lidar.2)))
  
  ## HexPts do not spatially match lidar plots
  ## function to allocate coordinates of x (sample points) those of the nearest y (lidar plots)
  ## requires data to be in ppp format 
  ## also limit assignment to lidar plots within 25m of sample points to remove sample points outside of forest dist and undist condition classes 
  x_ppp <- as.ppp(HexPts.250m.SPDF)  ## ppp object (point pattern) of the spatstat library
  y_ppp <- as.ppp(lidar.2.fc.dist.undist)
  nearest_y_from_x <- nncross(x_ppp, y_ppp)
  x_new_coords <- y_ppp[ nearest_y_from_x$which, ]  ## nearest_y_from_x$which contains the indices to select the closest points in y_ppp
  x_new <- as.SpatialPoints.ppp(x_new_coords) ## converting to an object of the class SpatialPoints
  FCID <- x_new_coords$marks$FCID ## get the IDs 
  x_new_dataframe <- data.frame(nearest_y_from_x, FCID)  ## and create a dataframe
  x_new <- SpatialPointsDataFrame(x_new, x_new_dataframe) ## add dataframe table to x_new
  x_new@proj4string <- lidar.2.fc.dist.undist@proj4string  ## assign same projection as lidar.2.fc.dist.undist
  HexPts.250m.Snap <- subset(x_new, x_new$dist < params0$min.snap.dist) ## limits assignment to lidar plots within 25m of sample points to remove sample points outside of forest dist and undist condition classes
  
  ## create polygons (used for filtering) to match 3 * 3 window of lidar plots
  radius <- 37.5  ## size of half the side of the sampling polygon unit (37.5 --> 75 m polygon in original Harold's work)
  xy.coords <- data.frame(unlist(coordinates(HexPts.250m.Snap)))
  yPlus <- xy.coords$my+radius
  xPlus <- xy.coords$mx+radius
  yMinus <- xy.coords$my-radius
  xMinus <- xy.coords$mx-radius
  FCID <- HexPts.250m.Snap@data$FCID
  square <- cbind(xMinus,yPlus, xPlus,yPlus, xPlus,yMinus, xMinus,yMinus, xMinus,yPlus)
  ## loop to create polygons from sample points
  a <- vector('list', length(2))
  for (i in 1:length(FCID)) {     
    a[[i]]<-Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), FCID[i]) 
  }
  polys.250m <- SpatialPolygons(a, proj4string=CRS(proj4string(HexPts.250m.Snap)))
  polys.250m.df <- SpatialPolygonsDataFrame(polys.250m, data.frame(id=FCID, row.names=FCID))

  ## keep only forested disturbed and undisted lidar plots within polygons
  lidar.2.fc.dist.undist.polys.250m <- lidar.2.fc.dist.undist[!is.na(over(lidar.2.fc.dist.undist, as(polys.250m.df, "SpatialPolygons"))),]
  
  ## get number of lidar plots in each polygon and keep only polygons with 9 plots 
  nptspoly.250m <- over(lidar.2.fc.dist.undist.polys.250m, polys.250m.df)
  table.250m <- data.frame(table(nptspoly.250m$id))
  freq.thresh <- 9  ## for this filtering step we work at the 3x3 (75m) spatial unit size, so 9 plots/polygon
  table.250m.subset <- subset(table.250m, Freq == freq.thresh)   
  
  ## keep only polygons with at least 9 lidar plots
  polys.250m.9pts.df <- polys.250m.df[polys.250m.df$id %in% table.250m.subset$Var1, ]

  ## keep lidar plots in sample polygons that have 9 plots
  lidar.2.fc.dist.undist.polys.250m.9pts <- lidar.2.fc.dist.undist[!is.na(over(lidar.2.fc.dist.undist,as(polys.250m.9pts.df,"SpatialPolygons"))),] 
  
  ## get polygon plot id for each lidar plot and rename it
  polyid <- over(lidar.2.fc.dist.undist.polys.250m.9pts, polys.250m.9pts.df)
  names(polyid) <- "POLY250ID"
  ## add polyid to data for lidar plots
  lidar.2.fc.dist.undist.polys.250m.9pts@data$POLY250ID <- polyid$POLY250ID
  
  ## subset to get centerpoint within each sample polygons
  lidar.2.fc.dist.undist.polys.250m.centerpts <- subset(lidar.2.fc.dist.undist.polys.250m.9pts, FCID==POLY250ID)
  ## export center lidar plots within sample polygons
  writeOGR(lidar.2.fc.dist.undist.polys.250m.centerpts, wkg_dir, paste(zone,"_pt_centerpt", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)

  ## build polygons used for sampling in prog2 
  if (paramsGL$unit.size == 75) {  ## if the chosen unit size is 75m the sampling polygon corresponds to the filtering polygon
    polys.250m.sampling.df <- polys.250m.9pts.df   
  } else if (paramsGL$unit.size == 25) {      ## else if it is 25m then starting from the center of 3x3 polygon build a 1x1 polygon covering just the central plot to be used to sample the raster layers in prog2
    radius <- 12.5
    xy.coords <- data.frame(unlist(coordinates(lidar.2.fc.dist.undist.polys.250m.centerpts)))
    colnames(xy.coords) <- c("mx", "my")
    yPlus <- xy.coords$my+radius
    xPlus <- xy.coords$mx+radius
    yMinus <- xy.coords$my-radius
    xMinus <- xy.coords$mx-radius
    FCID <- lidar.2.fc.dist.undist.polys.250m.centerpts@data$FCID
    square <- cbind(xMinus,yPlus, xPlus,yPlus, xPlus,yMinus, xMinus,yMinus, xMinus,yPlus)
    ## loop to create polygons from sample points
    a <- vector('list', length(2))
    for (i in 1:length(FCID)) {    
      a[[i]]<-Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), FCID[i]) 
    }
    polys.250m.sampling <- SpatialPolygons(a, proj4string=CRS(proj4string(polys.250m.9pts.df)))
    polys.250m.sampling.df <- SpatialPolygonsDataFrame(polys.250m.sampling, data.frame(id=FCID, row.names=FCID))
  }
  writeOGR(polys.250m.sampling.df, wkg_dir, paste(zone,"_poly_sampling",sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## export training polygons and lidar plots
  writeOGR(lidar.2.fc.dist.undist.polys.250m.9pts,  wkg_dir, paste(zone,"_pt_9plots_filtering",sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)  #  final layer with the 3x3 LiDAR plots selected 

#   nr.samples.per.zone[z] <- nrow(polys.250m.9pts.df@data)   ## to be used with sequential for
  nrow(polys.250m.9pts.df@data)
  
}

print(paste("Total nr. of samples:", sum(nr.samples.per.zone)))

stopCluster(cl)

## clock global time
toc <- proc.time()-tic[3]
print(paste("Prog0, total elapsed time:",seconds_to_period(toc[3])))
