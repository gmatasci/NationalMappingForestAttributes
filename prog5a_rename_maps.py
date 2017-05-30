
from __future__ import division  # to allow floating point divisions
import os
import sys
import glob
import logging
import shutil     # to remove folders and their contents
import gc       # to manually run garbage collection
import arcpy
import numpy as np
import pandas as pd
from collections import OrderedDict   # to have dictionaries whose order does not change

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*
from Functions_NatMapping_Python import*

print(python_info())

## PARAMETERS -------------------------------------------------------------------

params = {}

# params['zones'] = np.arange(7, 22+1, 1)
params['zones'] = np.array([10])

# params['subzones'] = ['S', 'N']
params['subzones'] = ['N']


params['renaming_dict'] = OrderedDict([('gross_stem_volume', 'gross_stem_volume_WRONG'), ('gross_stem_volume_NEW', 'gross_stem_volume'), ('total_biomass', 'total_biomass_WRONG'), ('total_biomass_NEW', 'total_biomass')])

root_dir = r'D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes'
ntems_dir = r'E:\NTEMS'

UTM_dir = os.path.join(ntems_dir, 'UTM_results')

# temp_dir = os.path.join(root_dir, 'wkg', 'temp')


## START ------------------------------------------------------------------------

print('Start')

start_time = tic()

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

log_file_path = os.path.join(root_dir, 'wkg_BOREAL', 'log_prog5a.txt')
if os.path.exists(log_file_path):
    os.remove(log_file_path)
logging.basicConfig(filename=log_file_path, level=logging.DEBUG)

# arcpy.env.workspace = temp_dir  # set temporary files directory (used when calling Reclassify(), etc.), could use arcpy.env.workspace = "in_memory"

## Delete if existing and create temp directory
# if os.path.exists(temp_dir):
#     shutil.rmtree(temp_dir)
# os.makedirs(temp_dir)
# arcpy.env.workspace = temp_dir  # set temporary files directory (used when calling Reclassify(), etc.)

for idx_zone, zone in enumerate(params['zones']):
    for idx_subzone, subzone in enumerate(params['subzones']):

        UTMz_name = '%d%s' % (zone, subzone)

        map_dir = os.path.join(UTM_dir, 'UTM_' + UTMz_name, 'YAI')

        for old_name, new_name in params['renaming_dict'].iteritems():

            # print('\t\t %s' % targ)

            try:
                file_key = os.path.join(map_dir, '*_' + old_name + '_YAI_2010.dat')
                map_path_old = glob.glob(file_key)[0]
                map_path_new = map_path_old.replace(old_name, new_name)
            except:
                logging.info('No such file: %s' % (file_key))
                continue

            arcpy.Rename_management(map_path_old, map_path_new)