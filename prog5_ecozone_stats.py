"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: prog5_ecozone_stats.py                            
Objective: Compute descriptive stats by ecozone over the entire boreal for all the mapped attributes
"""

## TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - reset all parameters to full run values

## SOLVED:
# -V Artifacts with biomass = 0 on UTM7S: issue happening also on other UTMzones. Has to do with the way topo layers were produced:
#   layers are larger by one pixel thus there is a one-pixel diagonal shift which results in zeros bc imputation has not even been asked to predict there
#   -- solved by Tx but we do not retrain the model, we just reapply it on the newly computed UTM zone layers
# -V substitute folders in E:\NTEMS\UTM_results with newly computed maps for UTM zones with extra pixel issues

## Failed attempt at making arcpy available to a Python version (3.5) coming with Anaconda: cannot import archook module
# sys.path.append(r'D:/Research/MyPythonModules/GitHubPython_modules')
# import archook #The module which locates arcgis
# archook.get_arcpy()
# import arcpy

## IMPORT MODULES ---------------------------------------------------------------

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

params['zones'] = np.arange(7, 22+1, 1)
# params['zones'] = np.array([17, 18])

params['subzones'] = ['S', 'N']
# params['subzones'] = ['S']

# params['ecozone_labels'] = OrderedDict([('Boreal Cordillera', 3), ('Boreal Plains', 4), ('Boreal Shield East', 5),
#                                        ('Boreal Shield West', 6), ('Hudson Plains', 7), ('Taiga Cordillera', 15),
#                                        ('Taiga Plains', 16), ('Taiga Shield East', 17), ('Taiga Shield West', 18)])
# params['ecozone_labels'] = OrderedDict([('Boreal Cordillera', 3), ('Boreal Plains', 4)])
params['ecozone_labels'] = OrderedDict([('Taiga Shield West', 18)])

params['forest_recoding'] = [[20, 0], [31, 0], [32, 0], [33, 0], [40, 0], [50, 0], [80, 0], [100, 0],
                             [81, 1], [210, 1], [220, 1], [230, 1]]   ## format for Arcpy's Reclassify
# params['forest_recoding_dict'] = {20: 0, 31: 0, 32: 0, 33: 0, 40: 0, 50: 0, 80: 0, 100: 0,
#                              81: 1, 210: 1, 220: 1, 230: 1}   ## dictionary format for reclassify_raster_numpy()

params['ecozone_NFI_stats'] = pd.DataFrame(columns=['biom_M_t', 'non_treed_biom_M_t', 'forest_K_ha'])
params['ecozone_NFI_stats'].loc['Atlantic Maritime'] = 1298.77, 6.39, 16295.67
params['ecozone_NFI_stats'].loc['Boreal Cordillera'] = 1354.87, 0.88, 19116.46
params['ecozone_NFI_stats'].loc['Boreal Plains'] = 2920.44, 14.86, 38454.65
params['ecozone_NFI_stats'].loc['Boreal Shield'] = 10041.23, 64.90, 131274.73
params['ecozone_NFI_stats'].loc['Hudson Plains'] = 221.97, 2.65, 9857.72
params['ecozone_NFI_stats'].loc['Taiga Cordillera'] = 494.22, 8.86, 6442.88
params['ecozone_NFI_stats'].loc['Taiga Plains'] = 2780.50, 0.16, 33601.17
params['ecozone_NFI_stats'].loc['Taiga Shield'] = 2246.74, 14.92, 46292.88

# params['targ_fact'] = OrderedDict([('elev_mean', 1000), ('elev_stddev', 1000), ('elev_cv', 1000), ('elev_p95', 1000),
#                                    ('percentage_first_returns_above_2m', 100), ('percentage_first_returns_above_mean', 100),
#                                    ('loreys_height', 1000), ('basal_area', 100), ('gross_stem_volume', 100), ('total_biomass', 10)])
params['targ_fact'] = OrderedDict([('elev_mean', 1000), ('elev_p95', 1000), ('percentage_first_returns_above_2m', 100),
                                   ('loreys_height', 1000), ('basal_area', 100), ('gross_stem_volume', 100), ('total_biomass', 10)])

aberr_zero_flag = 2**16-2  # flag with 65534 any aberrant attribute value (not to change the type of the array from uint16)
nr_landsat_pix_2_ha = 30 * 30 / 10000

root_dir = r'D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes'
ntems_dir = r'E:\NTEMS'

results_dir = os.path.join(root_dir, 'results', 'StatsByEcozone')
UTM_dir = os.path.join(ntems_dir, 'UTM_results')
ecozones_dir = os.path.join(UTM_dir, 'UTMrasters_Ecozones')
VLCE_dir = os.path.join(UTM_dir, 'UTMrasters_VLCE', '2010_overlap_free')
non_overl_dir = os.path.join(UTM_dir, 'UTMrasters_non_overlapping_areas')
boreal_dir = os.path.join(UTM_dir, 'UTMrasters_Brandt_Boreal')
temp_dir = os.path.join(root_dir, 'wkg', 'temp')


## START ------------------------------------------------------------------------

print('Start')

start_time = tic()

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

log_file_path = os.path.join(root_dir, 'wkg', 'log_prog4c.txt')
if os.path.exists(log_file_path):
    os.remove(log_file_path)
logging.basicConfig(filename=log_file_path, level=logging.DEBUG)

arcpy.env.workspace = temp_dir  # set temporary files directory (used when calling Reclassify(), etc.), could use arcpy.env.workspace = "in_memory"

## Delete if existing and create temp directory
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)
os.makedirs(temp_dir)
arcpy.env.workspace = temp_dir  # set temporary files directory (used when calling Reclassify(), etc.)

## Create results directory
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

## Formula computing average biomass in the ecozones, subtracting not-treed biomass from total above-ground biomass (million tonnes) on forest land
avg_biom_t_ha = ((params['ecozone_NFI_stats'][['biom_M_t']].as_matrix() - params['ecozone_NFI_stats'][['non_treed_biom_M_t']]).as_matrix() * 10**6) / (params['ecozone_NFI_stats'][['forest_K_ha']].as_matrix() * 10**3)
params['ecozone_NFI_stats']['avg_biom_t_ha'] = avg_biom_t_ha
params['ecozone_NFI_stats'].to_csv(os.path.join(results_dir, 'NFI_ecozone_stats.csv'))  # save in a separate csv by ecozones


## LOOP THROUGH ECOZONES ------------------------------------------------------------------------

ecoz_codes_in_UTMz_dict = {}
ecoz_stats_dict = {}
UTMz_aberr_zeros_dict = {}
first_ecoz = True
for ecoz_name, ecoz_code in params['ecozone_labels'].items():

    print('Ecozone: %s' % (ecoz_name))

    start_time_ecoz = tic()

    ecoz_samples_arr = np.ndarray(shape=(0,len(params['targ_fact'])), dtype=np.dtype('uint16'))
    tot_pix_ecoz = 0
    tot_nr_bor_forest_sampl_ecoz = 0
    for idx_zone, zone in enumerate(params['zones']):
        for idx_subzone, subzone in enumerate(params['subzones']):

            UTMz_name = '%d%s' % (zone, subzone)

## CHECK IF ECOZONE IN UTM ZONE ---------------------------------------------------------------

            ## Check if ecozone file exists to skip UTM zones for which we do not produce any maps
            try:
                file_key = os.path.join(ecozones_dir, 'UTM'+UTMz_name+'*.dat')
                ecozones_path = glob.glob(file_key)[0]
            except:
                if first_ecoz:
                    logging.info('No such file: %s' % (file_key))
                continue

            ## Iterations 2 and onwards: check right away if ecozone code is within the image
            if not first_ecoz:  # if it is not the 1st iteration ecoz_codes_in_UTMz_dict already exists
                if ecoz_code not in ecoz_codes_in_UTMz_dict[UTMz_name]:  # so we can check right away, without even reading the ecozone raster, if the current ecozone touches this UTM zone
                    continue

            ## Read ecozone file
            ecozones = arcpy.Raster(ecozones_path)
            ecozones_arr = arcpy.RasterToNumPyArray(ecozones)
            del ecozones

            ## Save unique ecozone codes found in the raster (done only for the first ecozone)
            if first_ecoz:
                unique_codes = np.unique(ecozones_arr)
                ecoz_codes_in_UTMz_dict[UTMz_name] = np.delete(unique_codes, np.where(unique_codes == 0))
                if ecoz_code not in ecoz_codes_in_UTMz_dict[UTMz_name]:  # if not found within the list of codes go to the next iteration
                    continue

## CREATE MASKS ------------------------------------------------------------------------

            print('\t UTM zone: %s' % (UTMz_name))

            ## Read non-overlapping mask
            try:
                file_key = os.path.join(non_overl_dir, 'Overlap_free_mask_UTM' + UTMz_name + '.dat')
                nonoverl_path = glob.glob(file_key)[0]
            except:
                if first_ecoz:
                    logging.info('No such file: %s' % (file_key))
                continue
            nonoverl = arcpy.Raster(nonoverl_path)
            mask_nonoverl_binary_arr, spatial_info_UTM = raster_2_array(nonoverl)  # store spatial infos only once, as layers are co-registered
            del nonoverl

            ## Create ecozone mask
            recoding_dict = {x: 0 for x in ecoz_codes_in_UTMz_dict[UTMz_name][ecoz_codes_in_UTMz_dict[UTMz_name] != ecoz_code]}
            recoding_dict[ecoz_code] = 1
            mask_ecoz_binary_arr = reclassify_raster_numpy(ecozones_arr, recoding_dict)
            del ecozones_arr
            mask_ecoz_nonoverl_binary_arr = mask_nonoverl_binary_arr * mask_ecoz_binary_arr
            del mask_nonoverl_binary_arr, mask_ecoz_binary_arr

            ## Read VLCE file
            try:
                file_key = os.path.join(VLCE_dir, '*_'+UTMz_name+'_*.dat')
                VLCE_path = glob.glob(file_key)[0]
            except:
                if first_ecoz:
                    logging.info('No such file: %s' % (file_key))
                continue
            VLCE = arcpy.Raster(VLCE_path)

            ## Create VLCE mask
            remap_VLCE = arcpy.sa.RemapValue(params['forest_recoding'])
            mask_VLCE_binary = arcpy.sa.Reclassify(VLCE, 'VALUE', remap_VLCE, 0)  # mask of 1s and 0s
            del VLCE
            mask_VLCE_binary_arr, _ = raster_2_array(mask_VLCE_binary)  # store spatial infos only once, as layers are co-registered
            del mask_VLCE_binary

            ## Read boreal mask
            try:
                file_key = os.path.join(boreal_dir, 'Boreal_UTM'+UTMz_name+'*.dat')
                boreal_path = glob.glob(file_key)[0]
            except:
                if first_ecoz:
                    logging.info('No such file: %s' % (file_key))
                continue
            mask_boreal_binary = arcpy.Raster(boreal_path)
            mask_boreal_binary_arr, _ = raster_2_array(mask_boreal_binary)

            ## Create complete mask
            ## With Numpy
            mask_complete_arr = mask_ecoz_nonoverl_binary_arr * mask_VLCE_binary_arr * mask_boreal_binary_arr
            del mask_VLCE_binary_arr, mask_boreal_binary_arr

            ## With ArcPy: not working for UTM zone 10S bc of weird error "Process finished with exit code -1073741819 (0xC0000005)"
            # mask_complete = mask_VLCE_binary * mask_ecoz_binary * mask_boreal_binary
            # mask_complete_arr = arcpy.RasterToNumPyArray(mask_complete)

            tot_pix_ecoz += np.count_nonzero(mask_ecoz_nonoverl_binary_arr)  # increment count of pixels falling in that ecozone
            del mask_ecoz_nonoverl_binary_arr
            nr_bor_forest_sampl_UTMz = np.count_nonzero(mask_complete_arr)  # save number of samples to initialize temp_arr
            tot_nr_bor_forest_sampl_ecoz += nr_bor_forest_sampl_UTMz  # increment count of non-masked pixels

## LOOP THROUGH ATTRIBUTES ------------------------------------------------------------------------

            map_dir = os.path.join(UTM_dir, 'UTM_'+UTMz_name, 'YAI')
            temp_arr = np.empty((nr_bor_forest_sampl_UTMz, len(params['targ_fact']),)).astype(np.dtype('uint16'))
            for idx_targ, targ in enumerate(params['targ_fact']):

                # print('\t\t %s' % targ)

                try:
                    file_key = os.path.join(map_dir, '*_'+targ+'_*.dat')
                    map_path = glob.glob(file_key)[0]
                except:
                    if first_ecoz:
                        logging.info('No such file: %s' % (file_key))
                    continue
                map = arcpy.Raster(map_path)

                map_arr, _ = raster_2_array(map)
                map_w_flags_arr = reclassify_raster_numpy(map_arr, {0: aberr_zero_flag})
                del map_arr
                map_masked_arr = map_w_flags_arr * mask_complete_arr
                del map_w_flags_arr

                dict_key = '%s_%s' % (UTMz_name, targ)
                if dict_key not in UTMz_aberr_zeros_dict:
                    UTMz_aberr_zeros_dict[dict_key] = np.count_nonzero(map_masked_arr==aberr_zero_flag)

                temp_arr[:, idx_targ] = map_masked_arr[map_masked_arr != 0]  # at this point the only zeroes in the map are masked pixels (no forest, no boreal)
                del map_masked_arr

            del mask_complete_arr,
            ecoz_samples_arr = np.vstack((ecoz_samples_arr, temp_arr))
            del temp_arr
            gc.collect()

## COMPUTE DESCRIPTIVE STATS ------------------------------------------------------------------------

    print('\t Computing stats...')

    ecoz_samples_df = pd.DataFrame(data=ecoz_samples_arr, columns=list(params['targ_fact'].keys()))
    del ecoz_samples_arr
    ecoz_stats = ecoz_samples_df[ecoz_samples_df!=aberr_zero_flag].describe(percentiles=[.1, .25, .5, .75, .9])
    ecoz_stats[1:] = ecoz_stats[1:]/params['targ_fact'].values()
    ecoz_stats.loc['HectaresBorealForest', :] = ecoz_stats.loc['count', params['targ_fact'].keys()[0]]*nr_landsat_pix_2_ha
    ecoz_stats.loc['HectaresEcozone', :] = tot_pix_ecoz*nr_landsat_pix_2_ha
    ecoz_stats.loc['TotNrBorealForestPixInEcoz', :] = tot_nr_bor_forest_sampl_ecoz  # should be equal to Count
    del ecoz_samples_df
    ecoz_stats_dict[ecoz_name] = ecoz_stats  # store descr. stats in a dictionary by ecozone
    ecoz_stats.to_csv(os.path.join(results_dir, '%s_stats.csv' % ecoz_name))  # save in a separate csv by ecozones

    logging.info('%s: %s' % (ecoz_name, toc(start_time_ecoz)))

    if first_ecoz:  # if we reached this line means first iteration was successful , so we change the flag to False
        first_ecoz = False

    gc.collect()

for key, value in sorted(UTMz_aberr_zeros_dict.items()):
    logging.info('%s: %d aberrant zero pixels' % (key, value))

print(toc(start_time))
