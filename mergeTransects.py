"""
Project Name: NationalImputationForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: mergeTransects.py                        
Objective: Merge into a single shp all the individual shps containing the transects in each UTM zone
"""

import arcpy.mapping, sys, os, glob, fnmatch

from arcpy import env
from arcpy.sa import *

matches = []

rootPath = 'E:\LOP_data\LOP_transects'

counter = 0
for _, _, files in os.walk(rootPath):
    for UTMzone in range(8,25):
        for latZone in ["S", "N"]:
            fileName = "UTM%d%s_trsct.shp" % (UTMzone, latZone)
            if fileName in files:
                matchIn = os.path.join(rootPath, fileName)
                matchOut = os.path.join(rootPath, "Reproj_NAD1983LCC", fileName)
                matches.append (matchOut)
                coordinateSystem = "PROJCS[\"NAD_1983_Lambert_Conformal_Conic\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-95.0],PARAMETER[\"Standard_Parallel_1\",77.0],PARAMETER[\"Standard_Parallel_2\",49.0],PARAMETER[\"Latitude_Of_Origin\",49.0],UNIT[\"Meter\",1.0]]"                      
                arcpy.Delete_management(matchOut)
                arcpy.Project_management(matchIn, matchOut, coordinateSystem)
                counter = counter + 1
    break  ## to prevent the loop from entering the 2nd round that goes into the subdirectory
                   
mergedFileOut = r"E:\LOP_data\LOP_transects\CAN_trsct_reproj_NAD1983LCC.shp"
arcpy.Delete_management(mergedFileOut)
arcpy.Merge_management(matches, mergedFileOut)