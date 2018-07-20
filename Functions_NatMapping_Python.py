"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: Functions_NatMapping_Python.py                          
Objective: Define Python functions to make available to all py scripts.
"""

import arcpy
import numpy as np

def raster_2_array(raster_obj):
    """
    Convert an arcpy raster object into a numpy array

    Parameters
    ----------
    raster_obj: raster object
        Arcpy raster object to convert

    Returns
    -------
    raster_arr: numpy array
        Array with raster values
    spatial_info: dictionary
        Dictionary with the spatial metadata needed to georeference an array as a raster
    """
    dsc = arcpy.Describe(raster_obj)
    sr = dsc.SpatialReference
    ll = arcpy.Point(dsc.Extent.XMin, dsc.Extent.YMin) # store lower left corner to ensure the saved raster is put the right place
    raster_arr = arcpy.RasterToNumPyArray(raster_obj)
    try:
        spatial_info = {'lower_left': ll, 'spatial_reference': sr,
                 'cell_width': dsc.meanCellWidth, 'cell_height': dsc.meanCellHeight}
    except:
        spatial_info = {'lower_left': ll, 'spatial_reference': sr}
    return raster_arr, spatial_info

def array_2_raster(raster_arr, spatial_info):
    """
    Convert a numpy array into an arcpy raster object

    Parameters
    -------
    raster_arr: numpy array
        Array with raster values
    spatial_info: dictionary
        Dictionary with the spatial metadata needed to georeference an array as a raster

    Returns
    -------
    raster_obj: raster object
        Converted arcpy raster object
    """
    raster_obj = arcpy.NumPyArrayToRaster(raster_arr, spatial_info['lower_left'],
                                          spatial_info['cell_width'], spatial_info['cell_height'])
    arcpy.DefineProjection_management(raster_obj, spatial_info['spatial_reference'])
    return raster_obj

def group_arcmap_layers(mxd_in, df_in, empty_gr_lyr_path, grouping_dict):
    """
    Group ArcMap layers

    Parameters
    -------
    mxd_in: string
        Path to existing MXD
    df_in: string
        Name of dataframe in MXD
    empty_gr_lyr_path: string
        Path to empty group layer file .lyr
    grouping_dict: dictionary
        Dictionary with pairs of group name and corresponding file key (for glob.glob()) of the layers to group together

    Returns
    -------
    message: string
    """
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

def reclassify_raster_numpy(raster_arr, recode_dict):
    reclass_raster_arr = np.copy(raster_arr)
    for old, new in recode_dict.iteritems():
        reclass_raster_arr[raster_arr == old] = new
    return reclass_raster_arr



