"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca)
File Name: exportTimeSeries_CAN_YAI_1984_2012.py                        
Objective: Export to PNG the time-series of some specified attributes of interest zoomed to a given bookmark 
"""

# TO DO

# - che ck code with the addition of the landsat image
# - add difference between 2012 and 1984 and export it


## IMPORT MODULES ---------------------------------------------------------------

import arcpy, os, sys, glob
import numpy as np
from itertools import compress

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*

## PARAMETERS  ------------------------------------------------

mxds_folder = r"D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR\MXDs"

root_dir = os.path.join(mxds_folder, "TimeSeries", "FireHarvest_CAN_timeseries_1984_2012")

out_folder = os.path.join(root_dir, r"PNGs")

# input_mxd = os.path.join(root_dir, r'empty.mxd')  ## path to empty mxd that has to have a scale bar (will adapt automatically) and bookmarks already set
input_mxd = os.path.join(root_dir, r'empty_nothing.mxd')  ## path to empty mxd that has to have a scale bar (will adapt automatically) and bookmarks already set

output_mxd = os.path.join(root_dir, r"FireHarvest_CAN_timeseries_1984_2012.mxd")  ## path to final mxd to save results

params = {}

params['build_MXD'] = True
# params['build_MXD'] = False

params['export_PNGs'] = True

# params['years'] = np.arange(1984, 2013, 1)
params['years'] = np.arange(1984, 1986, 1)

# params['targ_fact'] = {'elev_p95': "elev_p95", "cover_2m": 'percentage_first_returns_above_2m', "stem_volume": 'gross_stem_volume', "ag_biomass": 'total_biomass'}
params['targ_fact'] = {"stem_volume": 'gross_stem_volume', "ag_biomass": 'total_biomass'}

in_folder = r"E:\NTEMS\UTM_results\YAI_temporal_maps"

mask_file_base_name = "LC_Class_HMM_*.dat"  ## UTMzone in old filename by TX
landsat_image_base_name = "SRef_UTM*_1984_proxy_v2.dat"

# mapping_areas = {"11H": "11S", "11HH": "11S", "15F": "15S", "18H":"18S", "13F":"13S"}
mapping_areas = {"11H": "11S", "11HH": "11S", "13F":"13S"}

lyrs_dir = os.path.join(mxds_folder, "LYRs")
mask_lyr_path = os.path.join(lyrs_dir, r'Mask_4_TreedClasses.lyr')
landsat_lyr_path = os.path.join(lyrs_dir, r'Landsat_RGB543.lyr.lyr')
attr_lyr_srckey = os.path.join(lyrs_dir, '*_contin_*.lyr')

## START --------------------------------------------------------------------------

print('exportTimeSeries_CAN_YAI_1984_2012.py: started on %s' % time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime()))

start_time = tic()

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

## Get MXD and dataframe
mxd = arcpy.mapping.MapDocument(input_mxd)
df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]  ## assuming there is only 1 df (called Layer, the default)
mxd.activeView = df.name
mxd.title = df.name

## LOAD LAYERS AND APPLY SYMBOLOGY ----------------------------------------------

if params['build_MXD']:

    print('Loading layers and applying symbology...')

    for mapping_area, UTMzone in mapping_areas.iteritems():  ## iterate over mapping areas

        print('Area: '+mapping_area)

        in_dir_YAI = os.path.join(in_folder, mapping_area, "Results", "YAI_fitted")
        VLCE_dir = os.path.join(in_folder, mapping_area, "Results", "VLCE")
        landsat_dir = os.path.join(in_folder, mapping_area, "Results", "proxy_values")


        file_dirs = [in_dir_YAI, VLCE_dir, landsat_dir]

        for idx_dir, file_dir in enumerate(file_dirs):  ## iterate over the directories

            if file_dir == in_dir_YAI:
                files_to_add_path = []
                source_lyrs_dict = {}
                for targ, targ_long in params['targ_fact'].iteritems():
                    in_file_base_name_YAI = "UTM_%s_%s_YAI_*.dat" % (mapping_area, targ_long)
                    files_to_add_path.extend(glob.glob(os.path.join(file_dir, in_file_base_name_YAI)))
            elif file_dir == VLCE_dir:
                files_to_add_path = glob.glob(os.path.join(file_dir, mask_file_base_name))
            elif file_dir == landsat_dir:
                files_to_add_path = glob.glob(os.path.join(file_dir, landsat_image_base_name))


            for idx_file, file_to_add in enumerate(files_to_add_path):

                file = file_to_add

                if file_dir == VLCE_dir:  ## select forest mask layer file
                    if UTMzone in file_to_add:
                        file_path_new = file_to_add.replace(UTMzone, mapping_area)
                        arcpy.Rename_management(file_to_add, file_path_new)
                        file = file_path_new

                layer_to_add = arcpy.mapping.Layer(file)
                arcpy.mapping.AddLayer(df, layer_to_add, "TOP")

                if file_dir == VLCE_dir:  ## select forest mask layer file
                    source_lyr_path = mask_lyr_path
                elif file_dir == landsat_dir:
                    source_lyr_path = landsat_lyr_path
                elif file_dir == in_dir_YAI:  ## or select corresponding attribute layer file
                    logic_idx = [targ_long in layer_to_add.name for targ, targ_long in params['targ_fact'].iteritems()]
                    targ = list(compress(params['targ_fact'].values(), logic_idx))[0]
                    cand_lyr_path_list = glob.glob(attr_lyr_srckey)
                    source_lyr_path = [cand_lyr_path for cand_lyr_path in cand_lyr_path_list if targ in cand_lyr_path][0]

                source_lyr = arcpy.mapping.Layer(source_lyr_path)
                layer_to_update = arcpy.mapping.ListLayers(mxd, layer_to_add.name, df)[0]  ## redefine layer_to_update as an arcpy object
                arcpy.ApplySymbologyFromLayer_management(layer_to_update, source_lyr)  ## first apply the symbology with ApplySymbologyFromLayer_management (needs a path to the layer providing the symbology)
                arcpy.mapping.UpdateLayer(df, layer_to_update, source_lyr, True)  ## update symbology by using the one of source_lyr (has to be an arpy object and not just the path)

                arcpy.RefreshActiveView()
                arcpy.RefreshTOC()



    ## SAVE MXD -------------------------------------------------------

    print('Saving MXD...')

    mxd.saveACopy(output_mxd)

else:

    mxd = arcpy.mapping.MapDocument(output_mxd)
    df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]  ## assuming there is only 1 df (called Layer, the default)
    mxd.activeView = df.name
    mxd.title = df.name

## ZOOM ON BOOKMARKS AND EXPORT PNGs --------------------------------------

if params['export_PNGs']:

    print('Exporting PNGs...')

    ## Turn of all layers in list
    for layer in arcpy.mapping.ListLayers(mxd, '', df):
        layer.visible = False
    arcpy.RefreshActiveView()

    ## Loop over the bookmarks
    for bkmk in arcpy.mapping.ListBookmarks(mxd):

        if bkmk.name[0:3] not in mapping_areas.keys():
            continue

        print('Area: %s' % bkmk.name)

        ## Set a projection system for the dataframe matching the UTM zone of the bookmark/mapping area
        df.spatialReference = arcpy.SpatialReference("NAD 1983 UTM Zone %sN" % bkmk.name[0:2])
        df.extent = bkmk.extent  ## set extent of dataframe to extent of bookmark

        ## Create directory for output .png files
        png_dir = os.path.join(out_folder, bkmk.name)
        if not os.path.exists(png_dir):
            os.makedirs(png_dir)

        search_key = "SRef_UTM%s_1984_proxy_v2.dat" % (bkmk.name.split('_')[0])
        try:
            landsat_layer = arcpy.mapping.ListLayers(mxd, search_key, df)[0]
        except:
            continue
        landsat_layer.visible = True
        arcpy.RefreshActiveView()
        png_name_out = "%s_Landsat543_1984.png" % (bkmk.name)  # build dynamic filename with the equivalent of sprintf
        arcpy.mapping.ExportToPNG(mxd, os.path.join(png_dir, png_name_out))  ## export mxd (zoomed on bookmark extent, and with only one layer available) to png
        landsat_layer.visible = False

        ## Loop over year of interest
        for yr in params['years']:

            search_key = 'LC_Class_HMM_%s_%d.dat' % (bkmk.name.split('_')[0], yr)
            try:
                mask_layer = arcpy.mapping.ListLayers(mxd, search_key, df)[0]
            except:
                continue

            placement = "TOP"
            mask_layer.visible = True

            for targ, targ_long in params['targ_fact'].iteritems():

                search_key = "UTM_%s_%s_YAI_%d.dat" % (bkmk.name.split('_')[0], targ_long, yr)
                try:
                    map_layer = arcpy.mapping.ListLayers(mxd, search_key, df)[0]
                except:
                    continue

                placement = "BOTTOM"
                map_layer.visible = True
                arcpy.RefreshActiveView()

                png_name_out = "%s_%s_%d.png" % (bkmk.name, targ, yr)  # build dynamic filename with the equivalent of sprintf
                arcpy.mapping.ExportToPNG(mxd, os.path.join(png_dir, png_name_out))  ## export mxd (zoomed on bookmark extent, and with only one layer available) to png
                map_layer.visible = False

            mask_layer.visible = False


print('Total ' + toc(start_time))















