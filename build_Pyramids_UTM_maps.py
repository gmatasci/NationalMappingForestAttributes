## CODE INFOS -------------------------------------------------------------------

## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: build_Pyramids_UTM_maps.py                            
## Objective: Build pyramids for all forest attributes maps in each UTM zone saved in 'E:\NTEMS\UTM_results'

## IMPORT MODULES ---------------------------------------------------------------

import arcpy
import os
import sys
import glob
import logging
import multiprocessing
import numpy as np
from arcpy import env

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*

## PARAMETERS -------------------------------------------------------------------

# Set the current workspace
#env.workspace = r'E:\NTEMS\UTM_results\UTM_13S_percentage_first_returns_above_2m'
#env.workspace = r'E:\NTEMS\UTM_results\RF_corrected'
#env.workspace = r'E:\NTEMS\UTM_results\YAI_corrected'
root_dir = r'E:\NTEMS\UTM_results'

params = {}

params['zones'] = np.arange(7, 22+1, 1)
# params['zones'] = np.array([17, 18])
params['zones'] = params['zones'][params['zones'] != 11]

params['subzones'] = ['S', 'N']
# params['subzones'] = ['S']


## START ------------------------------------------------------------------------

print('Start')

start_time = tic()

log_file_path = os.path.join(root_dir, 'log_pyramids.txt')
if os.path.exists(log_file_path):
    os.remove(log_file_path)
logging.basicConfig(filename=log_file_path, level=logging.DEBUG)

for idx_zone, zone in enumerate(params['zones']):
    for idx_subzone, subzone in enumerate(params['subzones']):

        UTMz_name = '%d%s' % (zone, subzone)
        maps_dir = os.path.join(root_dir, 'UTM_' + UTMz_name, 'YAI')
        ## Check if ecozone file exists to skip UTM zones for which we do not produce any maps
        try:
            maps_dir = glob.glob(maps_dir)[0]
        except:
            logging.info('No such folder: %s' % (maps_dir))
            continue

        print('UTM zone: %s' % (UTMz_name))

        env.workspace = maps_dir
        dat_files_list = arcpy.ListFiles("*.dat")
        logging.info('Nr map files for %s: %d' % (UTMz_name, len(dat_files_list)))

        ## Build pyramids for all .dat files in folder
        for layer in dat_files_list:
            arcpy.BuildPyramids_management(os.path.join(env.workspace, layer), skip_existing="SKIP_EXISTING")  ## do not create pyramids if they already exist

print(toc(start_time))
	
## Test with parallel approach (not exactly working!)
#pool = multiprocessing.Pool(24, None, None, 1)
#
#pool.map(arcpy.BuildPyramids_management, arcpy.ListFiles())
#
## closing pool and waiting for each task to be finished
#pool.close()
#pool.join()