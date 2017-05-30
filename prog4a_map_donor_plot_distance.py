"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: prog4a_map_donor_plot_distance.py
Objective: Computes the distance from each pixel to its donor plot (in Lat/Long degree units) and also saves Lat and Long of donor plot
"""

## TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - reset all parameters to full run values

## SOLVED:
#

## IMPORT MODULES ---------------------------------------------------------------

from __future__ import division  # to allow floating point divisions
import os
import sys
import glob
import logging
import shutil     # to remove folders and their contents
import gc       # to manually run garbage collection
#import arcpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import pandas as pd
import multiprocessing as mp
from functools import partial
from collections import OrderedDict   # to have dictionaries whose order does not change
from pyproj import Proj, transform
from scipy.spatial.distance import cdist
from math import radians, cos, sin, asin, sqrt

sys.path.append(r'D:/Research/MyPythonModules')
sys.path.append(r'D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code')
from mybasemodule import*
from Functions_NatMapping_Python import*

## Returns the distance in Km between two long, lat coordinates (or arrays of coordinates) in decimal degrees
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance (Haversine formula) between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371            ## Radius of earth in kilometers. Use 3956 for miles
    return c * r

## Function to compute distance from donor plot to pixels
def distance_maps(arr_in, X_trn_in):
    mapID_in = arr_in[:,:,0]
    coord_arr_in = arr_in[:,:,1:3]
    dist_out = np.empty(shape=(mapID_in.shape[0], mapID_in.shape[1], 3)).astype(np.dtype('uint32'))

    for idx in np.arange(len(X_trn_in)):

        if idx % 5000 == 0:
            print(idx)

        long = X_trn_in.loc[idx, 'Long']
        lat = X_trn_in.loc[idx, 'Lat']
        bool_arr = (mapID_in == idx+1)  ## to have indices matching, as in maps they start at 1
        XY = np.vstack((coord_arr_in[:, :, 0][bool_arr], coord_arr_in[:, :, 1][bool_arr])).transpose()

        dist = haversine(XY[:,0], XY[:,1], long, lat)   ## for euclidean distance in lat/long degree units: np.squeeze(cdist(XY, np.array([[long, lat]]), 'euclidean'))

        dist_out[:, :, 0][bool_arr] = np.round(np.absolute(long) * 100).astype(np.dtype('uint32'))   ## Longitude in map is to be divided by a factor of 100
        dist_out[:, :, 1][bool_arr] = np.round(lat * 100).astype(np.dtype('uint32'))
        dist_out[:, :, 2][bool_arr] = np.round(dist * 10).astype(np.dtype('uint32'))

    return dist_out

## Given pixel indice pair (i, j) in image units, it returns the (X, Y) coordinate pair in spatial units
def XYcoords_pixel(i, j):
    coordX = rasXmin + ((0.5 + j) * rasMeanCellWidth)
    coordY = rasYmax - ((0.5 + i) * rasMeanCellHeight)
    return coordX, coordY


if __name__ == '__main__':      ## to prevent the parallel workers from running the whole script but just letting them read the function sent to them by pool.map/pool.apply

    ## PARAMETERS -------------------------------------------------------------------

    params = {}

    params['test_mode'] = False
    # params['test_mode'] = True
    params['test_mode_nr_rows_raster'] = 2000
    params['test_mode_nr_rows_Xtrn'] = 60505

    # params['zones'] = np.arange(7, 22+1, 1)
    params['zones'] = np.array([9])

    # params['subzones'] = ['S', 'N']
    params['subzones'] = ['N']

    params['nr_tiles'] = mp.cpu_count()
    # params['nr_tiles'] = 8


    root_dir = r'D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes'
    res_dir = os.path.join(root_dir, 'results_BOREAL', 'plot_dist_maps')
    wkg_dir = os.path.join(root_dir, 'wkg_BOREAL')
    ntems_dir = r'E:\NTEMS'
    UTM_dir = os.path.join(ntems_dir, 'UTM_results')
    temp_dir = os.path.join(wkg_dir, 'temp')

    X_trn_val_path = os.path.join(wkg_dir, 'X_trn_val.csv')

    ## START ------------------------------------------------------------------------

    print(python_info())

    print('Start')

    start_time = tic()

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True

    log_file_path = os.path.join(wkg_dir, 'log_prog4a.txt')
    if os.path.exists(log_file_path):
        os.remove(log_file_path)
    logging.basicConfig(filename=log_file_path, level=logging.DEBUG)

    arcpy.env.workspace = temp_dir  # set temporary files directory (used when calling Reclassify(), etc.), could use arcpy.env.workspace = "in_memory"

    ## Delete if existing and create temp directory
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)
    arcpy.env.workspace = temp_dir  # set temporary files directory (used when calling Reclassify(), etc.)

    ## Delete if existing and create temp directory
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    ## Read old and new training sets and merge them
    X_trn_val = pd.read_csv(X_trn_val_path)
    X_trn = X_trn_val[X_trn_val['TV'] == "TRAINING"]
    X_trn.reset_index(drop=True, inplace=True)   ## to have a continuous index after having removed some rows

    if params['test_mode']:
        X_trn = X_trn.loc[0:params['test_mode_nr_rows_Xtrn']-1,:]

    ## Loop of the UTM zones
    for idx_zone, zone in enumerate(params['zones']):
        for idx_subzone, subzone in enumerate(params['subzones']):

            UTMz_name = '%d%s' % (zone, subzone)
            print('\t UTM zone: %s' % (UTMz_name))

            map_dir = os.path.join(UTM_dir, 'UTM_' + UTMz_name, 'YAI')
            file_key = os.path.join(map_dir, '*ID_YAI_2010.dat')
            mapID_path = glob.glob(file_key)[0]

            ## Get properties of the input raster
            inRasterDesc = arcpy.Describe(mapID_path)

            ## Coordinates of the lower left corner
            rasXmin = inRasterDesc.Extent.XMin
            rasYmax = inRasterDesc.Extent.YMax

            ## Cell size, raster size
            rasMeanCellHeight = inRasterDesc.MeanCellHeight
            rasMeanCellWidth = inRasterDesc.MeanCellWidth
            
            if params['test_mode']:
                nr_rows = params['test_mode_nr_rows_raster']
                nr_cols = params['test_mode_nr_rows_raster']
            else:
                nr_rows = inRasterDesc.Height
                nr_cols = inRasterDesc.Width

            ## Execute the function XYcoords_pixel over each coordinate to get grid of long/lat values
            inRasterXY_grids = np.fromfunction(XYcoords_pixel, (nr_rows, nr_cols))  ## outputs tuple of X and Y coords
            inRasterXY = np.vstack((inRasterXY_grids[0].reshape(1,nr_rows*nr_cols), inRasterXY_grids[1].reshape(1,nr_rows*nr_cols))).transpose()  ## reshaped as an array with nr_pix rows and 2 cols to be input to transform()
            del inRasterXY_grids

            inProj = Proj(proj="utm", zone=zone, datum='WGS84')  ## all these are equivalent for utm zone 7N   Proj(init='epsg:%d' % (32600+zone)),  Proj('+proj = utm + zone = 7 + datum = WGS84 + units = m + no_defs')
            outProj = Proj(proj='latlong', datum='WGS84')  ## equivalent to Proj('+proj=longlat')
            long, lat = transform(inProj, outProj, inRasterXY[:,0], inRasterXY[:,1])
            del inRasterXY

            coord_arr = np.dstack((long.reshape(nr_rows, nr_cols), lat.reshape(nr_rows, nr_cols)))   ## 3D array with 1st layer the grid of longitude values and 2nd layer the grid of latitude values
            del long, lat

            mapID = arcpy.Raster(mapID_path)

            mapID_arr, sp_info = raster_2_array(mapID)

            if params['test_mode']:
                # mapID_arr = mapID_arr[5000:5000+nr_rows, 2000:2000+nr_cols]
                mapID_arr = mapID_arr[0:nr_rows, 0:nr_cols]
                coord_arr = coord_arr[0:nr_rows, 0:nr_cols, :]

            stack_arr = np.dstack((mapID_arr, coord_arr))
            del mapID_arr, coord_arr

            stack_arr_tiles = np.array_split(stack_arr, params['nr_tiles'])  ## splits into evenly-sized tiles (rectangular blocks)

            ## SEQUENTIAL: TO COMMENT
            # dist_arr = distance_maps(arr_in=stack_arr, X_trn_in=X_trn)
            ## SEQUENTIAL: TO COMMENT

            del stack_arr

            ## PARALLEL: TO UNCOMMENT
            nr_cores = min(len(stack_arr_tiles), mp.cpu_count())
            pool = mp.Pool(processes=nr_cores)             ## Create a multiprocessing Pool
            partial_distance_maps = partial(distance_maps, X_trn_in=X_trn)          ## partial() is used to specify which arguments remain constant (Y_trn)...
            dist_arr_tiles = pool.map(partial_distance_maps, stack_arr_tiles)      ## ...so that in pool.map() call only the list to iterate over/break into chunks is provided (map_stack_tiles)
            pool.close()            ## terminate worker processes
            pool.join()
            dist_arr = np.concatenate(dist_arr_tiles)        ## reput together the output tiles
            del dist_arr_tiles
            ## PARALLEL: TO UNCOMMENT

            long_map = array_2_raster(dist_arr[:,:,0], sp_info)      ## convert map back to a geospatial array with info saved in sp_info
            long_map_path = os.path.join(res_dir, 'Donor_plot_long_UTM_' + UTMz_name + '.dat')
            arcpy.CopyRaster_management(long_map, long_map_path, '32_BIT_UNSIGNED')     ## save it as a .dat file
            arcpy.DefineProjection_management(long_map_path, sp_info['spatial_reference'])
            del long_map

            lat_map = array_2_raster(dist_arr[:,:,1], sp_info)
            lat_map_path = os.path.join(res_dir, 'Donor_plot_lat_UTM_' + UTMz_name + '.dat')
            arcpy.CopyRaster_management(lat_map, lat_map_path, '32_BIT_UNSIGNED')
            arcpy.DefineProjection_management(lat_map_path, sp_info['spatial_reference'])
            del lat_map

            dist_map = array_2_raster(dist_arr[:,:,2], sp_info)
            dist_map_path = os.path.join(res_dir, 'Donor_plot_dist_UTM_' + UTMz_name + '.dat')
            arcpy.CopyRaster_management(dist_map, dist_map_path, '32_BIT_UNSIGNED')
            arcpy.DefineProjection_management(dist_map_path, sp_info['spatial_reference'])
            del dist_map

            gc.collect()

    print('Total '+toc(start_time))

