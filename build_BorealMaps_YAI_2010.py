"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: build_BorealMaps_YAI_2010.py                            
Objective: Build MXD with the 10 forest attributes maps for 2010 masked on the boreal
"""

## IMPORT MODULES ---------------------------------------------------------------

from __future__ import division  # to allow floating point divisions
import os
import fnmatch
import sys
import glob
import logging
import shutil     # to remove folders and their contents
import gc       # to manually run garbage collection
import arcpy
import numpy as np
import pandas as pd
from collections import OrderedDict
from itertools import compress

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*
from Functions_NatMapping_Python import*

print(python_info())

## PARAMETERS -------------------------------------------------------------------

mxd_dir = r'D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes\MXDs'
proj_dir = os.path.join(mxd_dir, 'Boreal_YAI_2010')
lyrs_dir = os.path.join(mxd_dir, "LYRs")
UTM_dir = r'E:\NTEMS\UTM_results'
VLCE_dir = os.path.join(UTM_dir, 'UTMrasters_VLCE', '2010_overlap_free')

UTM_subdir_key = "UTM_*" ## construct search key with variable parts replaced by *
# UTM_subdir_key = "UTM_2*" ## construct search key with variable parts replaced by *

log_file_path = os.path.join(proj_dir, 'log_build_BorealMaps_YAI_2010.txt')

# input_mxd = os.path.join(proj_dir, r'empty_Boreal_YAI_2010.mxd')  ## path to empty mxd that has to have a scale bar (will adapt automatically) and bookmarks already set
# output_mxd =  os.path.join(proj_dir, r'Boreal_YAI_2010.mxd')  ## path to final mxd to save results
input_mxd = os.path.join(proj_dir, r'empty_Boreal_YAI_2010_inset.mxd')  ## path to empty mxd that has to have a scale bar (will adapt automatically) and bookmarks already set
output_mxd =  os.path.join(proj_dir, r'Boreal_YAI_2010_inset_green.mxd')  ## path to final mxd to save results
# output_mxd =  os.path.join(proj_dir, r'Boreal_YAI_2010_inset_redgreen.mxd')  ## path to final mxd to save results
# output_mxd =  os.path.join(proj_dir, r'Boreal_YAI_2010_inset_brownblue.mxd')  ## path to final mxd to save results

new_empty_gr_lyr_path = os.path.join(lyrs_dir, r'NewGroupLayer.lyr')  ## path to empty layer file to be loaded to initialize each group
mask_lyr_path = os.path.join(lyrs_dir, r'Mask_4_TreedClasses.lyr')
attr_lyr_srckey = os.path.join(lyrs_dir, '*_contin_*_green.lyr')
# attr_lyr_diverging_srckey = os.path.join(lyrs_dir, '*_contin_*_redgreen.lyr')
# attr_lyr_diverging_srckey = os.path.join(lyrs_dir, '*_contin_*_brownblue.lyr')

layer_ref_name = 'Inverse_Boreal_CAN'   ## name of shp with which to mask out the raster layers with the maps

params = {}

params['years'] = np.array([2010])

params['targ_names'] = ['elev_mean', 'elev_sd', 'elev_cv', 'elev_p95',
                        'cover_2m', 'cover_mean',
                        'loreys_height', 'basal_area', 'stem_volume', 'ag_biomass']
params['targ_src_keys'] = ['elev_mean', 'elev_stddev', 'elev_cv', 'elev_p95',
                           'percentage_first_returns_above_2m', 'percentage_first_returns_above_mean',
                           'loreys_height', 'basal_area', 'gross_stem_volume', 'total_biomass']

params['mask_dict'] = {'Forest mask': ['Nonm_overlap_LC_Class_HMM_', 'TOP']}

params['method'] = 'YAI'


## START ------------------------------------------------------

print('Start')

start_time = tic()

if os.path.exists(log_file_path):
    os.remove(log_file_path)
logging.basicConfig(filename=log_file_path, level=logging.DEBUG)

positions = ['BOTTOM'] * len(params['targ_names'])
grouping_dict = {n:[params['targ_src_keys'][i], positions[i]] for (i, n) in enumerate(params['targ_names'])}
grouping_dict.update(params['mask_dict'])

file_dirs = [os.path.join(dir, params['method']) for dir in glob.glob(os.path.join(UTM_dir, UTM_subdir_key))]

file_dirs.append(VLCE_dir)

## Get MXD and dataframe
mxd = arcpy.mapping.MapDocument(input_mxd)
df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]  ## assuming there is only 1 df (called Layer, the default)
mxd.activeView = df.name
mxd.title = df.name

bkmk = arcpy.mapping.ListBookmarks(mxd, "", df)[0]
df.extent = bkmk.extent   ## set extent of dataframe to extent of bookmark

## LOAD LAYERS AND APPLY SYMBOLOGY --------------------------------------

print('Loading layers and applying symbology...')

for idx_dir, file_dir in enumerate(file_dirs):  ## iterate over the directories

    arcpy.env.workspace = file_dir  ## set this directory as workspace on which arcpy.ListFiles() will be called on
    dat_files_list = arcpy.ListFiles('*.dat')
    logging.info('Nr files for %s: %d' % (file_dir, len(dat_files_list)))

    for layer in dat_files_list:

        layer_path = os.path.join(arcpy.env.workspace, layer)
        try:
            layer_to_add = arcpy.mapping.Layer(layer_path)
        except:
            logging.info('No such file: %s' % (layer_path))
            continue
        arcpy.mapping.AddLayer(df, layer_to_add, "TOP")

        if file_dir == VLCE_dir:   ## select forest mask layer file
            source_lyr_path = mask_lyr_path
        else:  ## or select corresponding attribute layer file
            logic_idx = [name in layer_to_add.name for i, name in enumerate(params['targ_src_keys'])]
            targ = list(compress(params['targ_src_keys'], logic_idx))[0]
            cand_lyr_path_list = glob.glob(attr_lyr_srckey)
            source_lyr_path = [cand_lyr_path for cand_lyr_path in cand_lyr_path_list if targ in cand_lyr_path][0]

        source_lyr = arcpy.mapping.Layer(source_lyr_path)
        layer_to_update = arcpy.mapping.ListLayers(mxd, layer_to_add.name, df)[0]  ## redefine layer_to_update as an arcpy object
        arcpy.ApplySymbologyFromLayer_management(layer_to_update, source_lyr)  ## first apply the symbology with ApplySymbologyFromLayer_management (needs a path to the layer providing the symbology)
        arcpy.mapping.UpdateLayer(df, layer_to_update, source_lyr, True)  ## update symbology by using the one of source_lyr (has to be an arpy object and not just the path)


## GROUPING AND REORDERING LAYERS ---------------------------------------------------------------

print('Grouping layers...')

message = group_arcmap_layers(mxd, df, new_empty_gr_lyr_path, grouping_dict)
print(message)

print('Reordering layers...')
layer_ref = arcpy.mapping.ListLayers(mxd, layer_ref_name, df)[0]
layer_to_move = arcpy.mapping.ListLayers(mxd, params['mask_dict'].keys()[0], df)[0]
arcpy.mapping.MoveLayer(df, layer_ref, layer_to_move, "After")


## SAVE MXD -------------------------------------------------------

print('Saving MXD...')

arcpy.RefreshActiveView()
arcpy.RefreshTOC()

mxd.saveACopy(output_mxd)

print('Total '+toc(start_time))


