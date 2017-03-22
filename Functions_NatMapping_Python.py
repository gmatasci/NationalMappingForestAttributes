"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: Functions_NatMapping_Python.py                          
Objective: Define Python functions to make available to all py scripts.
"""

import arcpy
import numpy as np

def raster_2_array(raster_obj):
    dsc = arcpy.Describe(raster_obj)
    sr = dsc.SpatialReference
    ll = arcpy.Point(dsc.Extent.XMin, dsc.Extent.YMin) # store lower left corner to ensure the saved raster is put the right place
    raster_arr = arcpy.RasterToNumPyArray(raster_obj)
    spatial_info = {'lower_left': ll, 'spatial_reference': sr,
                 'cell_width': dsc.meanCellWidth, 'cell_height': dsc.meanCellHeight}
    return raster_arr, spatial_info

def array_2_raster(raster_arr, spatial_info):
    raster_obj = arcpy.NumPyArrayToRaster(raster_arr, spatial_info['lower_left'],
                                          spatial_info['cell_width'], spatial_info['cell_height'])
    arcpy.DefineProjection_management(raster_obj, spatial_info['spatial_reference'])
    return raster_obj

def reclassify_raster_numpy(raster_arr, recode_dict):
    reclass_raster_arr = np.copy(raster_arr)
    for old, new in recode_dict.iteritems():
        reclass_raster_arr[raster_arr == old] = new
    return reclass_raster_arr

def group_arcmap_layers(mxd_in, df_in, empty_gr_lyr_path, grouping_dict):
    i = 0
    for name, srckey_placement_list in grouping_dict.items():
        srckey = srckey_placement_list[0]
        placement = srckey_placement_list[1]
        group_layer = arcpy.mapping.Layer(empty_gr_lyr_path)  ## add empty group layer created beforehand
        group_layer.name = name
        arcpy.mapping.AddLayer(df_in, group_layer, placement)  ## add group layer to dataframe with a given placement ('TOP', 'BOTTOM')
        group_layer = arcpy.mapping.ListLayers(mxd_in, name, df_in)[0]  ## redefine group_layer as being the group layer named by targ (to avoid error later on)
        list_layers_targ = arcpy.mapping.ListLayers(mxd_in, '*' + srckey + '*', df_in)  ## list layers starting with the name of the fire
        for layer in list_layers_targ:
            if not layer.isGroupLayer:  ## avoids taking group layers
                arcpy.mapping.AddLayerToGroup(df_in, group_layer, layer)
                arcpy.mapping.RemoveLayer(df_in, layer)  ## remove layer to avoid duplicates
        i += 1
    return 'Layers grouped!'