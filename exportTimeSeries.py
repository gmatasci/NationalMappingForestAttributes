import arcpy
import os
import numpy as np

#%% PARAMETERS  ------------------------------------------------

mxdPath = r"D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes\MXDs"
#mxdName = r"Elev_p95_timeseries_to1996_NEW.mxd"
mxdName = r"Elev_p95_timeseries_OLD.mxd"
#targ = "Cover_2m"
targ = "Elev_p95"
fileBaseName = "UTM_13S_elev_p95_RF_"
outFolder = "TimeSeries"
years = np.arange(1984, 2012, 1)


#%% START ------------------------------------------------

## Create directort for output .png files
outDirectory = os.path.join(mxdPath, outFolder)
if not os.path.exists(outDirectory):
    os.makedirs(outDirectory)

## Get MXD and dataframe
mxd = arcpy.mapping.MapDocument(os.path.join(mxdPath, mxdName))
df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]  ## assuming there is only 1 df (called Layer, the default)

## Turn of all lyrs in list
for lyr in arcpy.mapping.ListLayers(mxd, '', df):
    lyr.visible = False
arcpy.RefreshActiveView()

## Loop over the bookmarks
for bkmk in arcpy.mapping.ListBookmarks(mxd):  
    df.extent = bkmk.extent   ## set extent of dataframe to extent of bookmark
    ## Loop over the layers in df
    for lyr in arcpy.mapping.ListLayers(mxd, '', df): 
        if lyr.name[:-8] == fileBaseName:  ## remove last 8 characters of string ("1995.dat") and check if it matches the basename 
            yr = lyr.name[len(lyr.name)-8:len(lyr.name)-4]    ## save year 
            if any(int(yr) == years):  ## if year is within list of years of interest
                lyr.visible = True   ## set layer visible
                arcpy.RefreshActiveView()
                fileNameOut = "%s_%s_%s.png" % (targ, bkmk.name, yr)   # build dynamic filename with the equivalent of sprintf 
                arcpy.mapping.ExportToPNG(mxd, os.path.join(outDirectory, fileNameOut))  ## export mxd (zoomed on bookmark extent, and with only one layer available) to png
                lyr.visible = False  ## deactivate layer





