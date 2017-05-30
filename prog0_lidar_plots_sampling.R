## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog0_lidar_plots_sampling.R                          
## Objective: Define sampling grid over LiDAR transect

#### TO DO -------------------------------------------------------------------

## STILL TO DO:

# Prior to actual run:


## SOLVED:
# -V change parts with XXXXXXXXX TO CHANGE XXXXXXXXXXXXXXXXXX tag to add filtering wrt for_sum, etc.
 


#### READS/WRITES ------------------------------------------------------------

## READS:
# - "<UTMzone>_points.shp": point shp with LiDAR 25x25 m plots
# - "<UTMzone>_trsct.shp": polyline shp with center of transect line

## WRITES:
# - "<UTMzone>_lop_wkg_fc.shp": point shp with plots surrounded by forest in a 3x3 neighborhood
# - "<UTMzone>_lop_wkg_fc_dist_undist.shp": point shp with plots surrounded by forest and all either disturbed or undisturbed in a 3x3 neighborhood (MAIN OUTPUT WHOSE FCIDs REFER TO) 
# - "<UTMzone>_pt_centerpt.shp": (for prog1) point shp with center of hexagons used for sampling, after subsetting wrt forest/dist/undist (MAIN OUTPUT USED TO PRODUCE SHP "<UTMzone>_cpt_poly_250m_training_validation" IN THE NEXT STEP)
# - "<UTMzone>_poly_sampling.shp": (for prog1) polygon shp with sampling polygons (3x3 or 1x1) covering plots' squares, after subsetting wrt forest/dist/undist
# - "<UTMzone>_pt_9plots_filtering.shp": (for prog1) point shp with the LiDAR plots (3x3 = 9 plots) used for filtering

#### INIT --------------------------------------------------------------------

print('Prog0: defining sampling grid')

rm(list=ls()) ## clear all variables

param_file = "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)


#### SCRIPT SPECIFIC PARAMETERS -------------------------------------------

params0 <- list()

# params0$lidar.sources <- c("NONBOREAL")  ## options: c("BOREAL", "NONBOREAL") but dont use BOREAL yet, not to overwrite files used for 2010 mapping (1st study)
params0$lidar.sources <- c("BOREAL", "NONBOREAL")

params0$transect.names <- c("AFRF", "CR04", "CR08", "MK10", "Quesnel", "Tofino")
# params0$transect.names <- c("CR04", "CR08")

params0$dist.center.trans <- 300   ## threshold on distance from center of transect to eliminate plots whose acquisition was too skewed
params0$hexag.cell.size <- 250   ## spacing of the hexagonal sampling cells
# params0$hexag.cell.size <- 800   ## spacing of the hexagonal sampling cells

params0$min.snap.dist <- 25  ## limit distance used in function to allocate coordinates of sample points those of the nearest lidar plots

param_file_prog0 = file.path(base_wkg_dir, 'AllUTMzones_params0.Rdata', fsep = .Platform$file.sep) 
save(params0, file = param_file_prog0)

#### LOAD PACKAGES ----------------------------------------------------------

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


#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time

nr.samples.per.zone <- array(NaN, length(paramsGL$zones))  ## info to report (output of foreach)

## read Ecozone polyg shp
ecoz_dir <- file.path(data_dir, 'NTEMS', 'Cartography', 'Ecozones', fsep = .Platform$file.sep)
ecozones <- readOGR(dsn = ecoz_dir, layer='Ecozones_CBM_LCC')  ## do not load "_Clip1" bc then some parts of the transect will be outside of the polygons 

nr.samples.per.zone.full <- vector(mode="numeric", length=0)
for (ls in 1:length(params0$lidar.sources)) {
  
  lidar.source <- params0$lidar.sources[ls]
  
  if (lidar.source == "BOREAL") {
    zones <- paramsGL$zones
    pnt_dir <- file.path(LOP_dir,'LOP_attributed', fsep = .Platform$file.sep)
    nr.clusters <- length(paramsGL$zones)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end
  } else {
    zones <- params0$transect.names
    pnt_dir <- file.path(LOP_dir,'LOP_attributed_BC', fsep = .Platform$file.sep)
    nr.clusters <- length(params0$transect.names)  ## for parallel just uncomment the foreach line, the preceding lines and the stopCluster(cl) line at the end
  }

  cl <- makeCluster(nr.clusters)
  registerDoParallel(cl)
  ## when using foreach the last unassigned variable can be stored in a growing dataframe (the "nr.samples.per.zone" to which the foreach result is assigned here)
  nr.samples.per.zone <- foreach (z=1:length(zones), .combine='rbind', .packages=list.of.packages) %dopar% {   ## add .verbose=TRUE for more info when debugging
    ## WHEN NOT USING FOREACH, TESTING PHASE
  # for (z in 1:length(zones)) {   ## to be used when debugging (using foreach breakpoints don't work and line numbers of errors is not shown)
    ## WHEN NOT USING FOREACH, TESTING PHASE
    
    zone <- zones[z]
    wkg_dir <- file.path(base_wkg_dir, zone, fsep = .Platform$file.sep)
  
    if (! file.exists(wkg_dir)){dir.create(wkg_dir, showWarnings = F, recursive = T)}
  
#### READ & CHECK FILES -----------------------------------------------------
    
    ## read lidar plot points shp
    lidar.raw <-  readOGR(dsn = pnt_dir, layer = paste(zone,"_points", sep='')) ## "lidar" is of type SpatialPointsDataFrame that is an object of the S4 object system, readOGR belongs to rgdal: dsn = data source name (folder), layer = filename without extension
    
    ## filter duplicate coordinates
    duplicIndic <- duplicated(lidar.raw@coords, incomparables = FALSE)
    lidar <- lidar.raw[!duplicIndic, ]  ## lidar now contains only one point at each coordinate pair
    rm(lidar.raw)
    rm(duplicIndic)
    
    ## reproject in this UTM zone the Ecozone polyg shp
    ecozones.UTMreproj <- spTransform(ecozones, CRS(proj4string(lidar)))
    
    ## read the lidar transect polyline shp
    if (lidar.source == "BOREAL") {
      trsct_dir <- file.path(LOP_dir,'LOP_transects', fsep = .Platform$file.sep) 
      lidar.transect <- readOGR(dsn = trsct_dir, layer = paste(zone, "_trsct",sep=''))
    
      ## check if coordinate systems are the same between lidar plots point shapefile and transect line shapefile
      if ( !all(sapply(list(proj4string(lidar.transect), proj4string(ecozones.UTMreproj)), FUN = identical, proj4string(lidar))) ) {
        stop(sprintf("%s: projections of LiDAR points shapefile, LiDAR transect shapefile and reprojected Ecozone shapefile loaded in prog0 do not match.", zone))
      }
      
    } else {
      
      ## check if coordinate systems are the same between lidar plots point shapefile and transect line shapefile
      if ( !all(sapply(list(proj4string(ecozones.UTMreproj)), FUN = identical, proj4string(lidar))) ) {
        stop(sprintf("%s: projections of LiDAR points shapefile and reprojected Ecozone shapefile loaded in prog0 do not match.", zone))
      }
      
    }
    
    ## update Ecozone field based on Ecozone shapefile
    lidar@data$ECOZONE <- over(lidar, ecozones.UTMreproj)[, "ECOZONE_NA"]  ## over() outputs a df with the values of the fields in ecozones.UTMreproj at the locations of lidar
  
    
#### FILTER WRT NEAR_DIST, FOR_SUM & DIST_SUM  -----------------------------------
    
    if (lidar.source == "BOREAL") {
      ## subset based on distance to edge of lidar transect: keep points that are no farther than params0$dist.center.trans meters (set to 300m) from the center of transect
      lidar.2 <- subset(lidar, NEAR_DIST <= params0$dist.center.trans)  ## NEAR_DIST is an attribute of the layer indicating the distance from the center of the transect
      names(lidar.2)[which(names(lidar.2) == "Unique_ID")] <- "unique_id"
    } else {  ## no subsetting wrt angle for NONBOREAL acquisitions as they are not a transect but multiple overpasses
      lidar.2 <- lidar
    }
    
    ## subset to keep only 3x3 plot configurations in forested areas (for_sum == 9)
    lidar.2.fc <- subset(lidar.2, for_sum == 9) 
    writeOGR(lidar.2.fc , wkg_dir, paste(zone,"_lop_wkg_fc", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)  ## save it for future checks
    
    ## further subset to keep only 3x3 plot configurations in completely disturbed (dist_sum == 9) or undisturbed (dist_sum == 0) areas
    lidar.2.fc.dist.undist <- subset(lidar.2.fc, dist_sum == 9 | dist_sum == 0 )
    
    ## check for duplicates in the unique_id field
    if (!all(duplicated(lidar.2.fc.dist.undist@data$unique_id) == F)) {
      stop(sprintf("%s: unique_id field has duplicates", zone))
    }
      
    ## use field unique_id to avoid duplicates across UTM zones when populating the new reference field FCID
    lidar.2.fc.dist.undist@data$FCID <- lidar.2.fc.dist.undist@data$unique_id   ## add unique identifier FCID for each of the points selected
    writeOGR(lidar.2.fc.dist.undist, wkg_dir, paste(zone,"_lop_wkg_fc_dist_undist", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)   ## shapefile with points respecting all the conditions (totally disturbed or totally undisturbed forest regions)
    

  #### DEFINE AND SNAP HEXAGONAL GRID ----------------------------------------------
    
    ## subset coordinates by defined minimum distance between plots
    ## sample points along hexagon lattice with 250m spacing
    ## restrict sampling to within 400 m buffer of either side of lidar transect center 
    
    if (lidar.source == "BOREAL") {
      transect.buffer <- gBuffer(lidar.transect, width=400)
    } else {
      transect.buffer <- gBuffer(lidar.2.fc.dist.undist, width=50)  ## buffer 50 m around scan areas for NONBOREAL acquisitions
    }
    
    set.seed(paramsGL$global.seed)  ## spssample has a random component in where it places the grid, so needs to be initialized the same each time
    HexPts.250m <- spsample(transect.buffer, type = "hexagonal", cellsize = params0$hexag.cell.size)    ## cell size is set to 250m, that is Harold's final value
  
    ## generate hexagon unique IDs
    Hex250ID <- as.data.frame(seq(1, length(HexPts.250m), by=1))
    names(Hex250ID) <- "Hex250ID"  ## assign name to the column of the data frame
    HexPts.250m.SPDF <- SpatialPointsDataFrame(coordinates(HexPts.250m), data=Hex250ID, proj4string=CRS(proj4string(lidar)))
    
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
  
    
  #### CREATE FILTERING POLYGONS ------------------------------------------
    
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
      a[[i]]<- Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), FCID[i]) 
    }
    polys.250m <- SpatialPolygons(a, proj4string=CRS(proj4string(HexPts.250m.Snap)))
    polys.250m.df <- SpatialPolygonsDataFrame(polys.250m, data.frame(id=FCID, row.names=FCID))
  
    
  #### FILTER WRT 9PTS/POLYG -------------------------------------------------------
    
    ## keep only forested disturbed and undisturbed lidar plots within polygons
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
  
  #### GET CENTERPOINTS OF 9PTS AFTER FILTERING ---------------------------------------
    
    ## get polygon plot id for each lidar plot and rename it
    polyid <- over(lidar.2.fc.dist.undist.polys.250m.9pts, polys.250m.9pts.df)
    names(polyid) <- "POLY250ID"
    ## add polyid to data for lidar plots
    lidar.2.fc.dist.undist.polys.250m.9pts@data$POLY250ID <- polyid$POLY250ID
    
    ## subset to get centerpoint within each sample polygons
    lidar.2.fc.dist.undist.polys.250m.centerpts <- subset(lidar.2.fc.dist.undist.polys.250m.9pts, FCID==POLY250ID)
    ## export center lidar plots within sample polygons
    writeOGR(lidar.2.fc.dist.undist.polys.250m.centerpts, wkg_dir, paste(zone,"_pt_centerpt", sep=''), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
    
  #### CREATE SAMPLING POLYGONS BASED ON CENTERPOINTS  ---------------------------------
    
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
    
    nrow(polys.250m.9pts.df@data)
    
  ## WHEN NOT USING FOREACH, TESTING PHASE
  #   nr.samples.per.zone[z] <- nrow(polys.250m.9pts.df@data)
  ## WHEN NOT USING FOREACH, TESTING PHASE
    
  }

  stopCluster(cl)
  
  nr.samples.per.zone.full <- rbind(nr.samples.per.zone.full, nr.samples.per.zone)
  
}


#### PRINT LOGS ---------------------------------------------------------

print(paste("Total nr. of samples:", sum(nr.samples.per.zone.full)))

## clock global time
toc <- proc.time()-tic[3]
print(paste("Prog0, total elapsed time:",seconds_to_period(toc[3])))