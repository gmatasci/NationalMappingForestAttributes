"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: exportMaps_YAI_2010.py                            
Objective: Export to PNG the 2010 maps of some specified attributes of interest zoomed to a given bookmark 
"""

## IMPORT MODULES ---------------------------------------------------------------

import arcpy, os, sys, glob
import numpy as np

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*

#%% PARAMETERS  ------------------------------------------------

mxds_folder = r"D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes\MXDs"

root_dir = os.path.join(mxds_folder, "Boreal_YAI_2010")

input_mxd = os.path.join(root_dir, r'empty_Boreal_YAI_2010.mxd')  ## path to empty mxd that has to have a scale bar (will adapt automatically) and bookmarks already set
output_mxd =  os.path.join(root_dir, r"Boreal_YAI_2010.mxd")  ## path to final mxd to save results

params = {}

params['years'] = np.array([2010])

params['targs'] = ["elev_p95", "cover_2m", "ag_biomass"]
params['targ_names_lg'] = ["elev_p95", "percentage_first_returns_above_2m", "total_biomass"]

params['symbol_type'] = "continuous"
# params['symbol_type'] = "classified"

# params['class_break_values'] = [2000, 7000, 10000, 13000, 16000, 60000]
# params['class_break_labels'] = ["2 to 7 m", "7 to 10 m", "10 to 13 m", "13 to 16 m", "> 16 m"]

params['method'] = "YAI"

#params['export_pngs'] = False
params['export_pngs'] = True

in_folder = r"E:\NTEMS\UTM_results\UTM_13S"

in_folder_YAI = os.path.join(in_folder, r'YAI') 

mask_folder = r"E:\NTEMS\UTM_results\VLCE\13S_HMM"
mask_file_base_name = "LC_Class_HMM_13S"

lyrs_folder = os.path.join(mxds_folder, "LYRs")

elev_p95_lyr_contin_path = os.path.join(lyrs_folder, r'Elev_p95_contin_5to19m.lyr')  ## path to layer file providing the symbology for the RGB composite
elev_p95_lyr_class_path = os.path.join(lyrs_folder, r'Elev_p95_class_0to23m.lyr')  ## path to layer file providing the symbology for the RGB composite
cover_2m_lyr_contin_path = os.path.join(lyrs_folder, r'Cover_2m_contin_0to100pct.lyr')  ## path to layer file providing the symbology for the RGB composite
ag_biomass_lyr_contin_path = os.path.join(lyrs_folder, r'Ag_biomass_contin_3to300tha.lyr')  ## path to layer file providing the symbology for the RGB composite

mask_lyr_path = os.path.join(lyrs_folder, r'Mask_4_TreedClasses.lyr')  ## path to layer file providing the symbology for the categorical rasters (unique values)

out_folder = os.path.join(root_dir, r"pngs")


## START ------------------------------------------------------------------------

print('Start')

start_time = tic()

file_dirs = [in_folder_YAI, mask_folder]

## Create directory for output .png files
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

## Get MXD and dataframe
mxd = arcpy.mapping.MapDocument(input_mxd)
df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]  ## assuming there is only 1 df (called Layer, the default)
mxd.activeView = df.name
mxd.title = df.name

#%% LOAD LAYERS AND APPLY SYMBOLOGY ----------------------------

print('Loading layers and applying symbology...')

# if params['symbol_type'] == "continuous":
#     elev_p95_lyr_path = elev_p95_lyr_contin_path
# elif params['symbol_type'] == "classified":
#     elev_p95_lyr_path = elev_p95_lyr_class_path

years_to_add = params['years']

for idx_dir, file_dir in enumerate(file_dirs):  ## iterate over the directories
    
    arcpy.env.workspace = file_dir   ## set this directory as workspace on which arcpy.ListFiles() will be called on
    os.chdir(file_dir)

    for yr in years_to_add:

        if file_dir == in_folder_YAI:

            for idx_targ, targ in enumerate(params['targs']):

                if targ == "elev_p95":
                    source_lyr_path = elev_p95_lyr_contin_path
                elif targ == "cover_2m":
                    source_lyr_path = cover_2m_lyr_contin_path
                elif targ == "ag_biomass":
                    source_lyr_path = ag_biomass_lyr_contin_path

                source_lyr = arcpy.mapping.Layer(source_lyr_path)
                search_key = "UTM_13S_%s_YAI_%d.dat" % (params['targ_names_lg'][idx_targ], yr)  ## construct search key with variable parts replaced by *

                file_to_add_path = glob.glob(search_key)[0]  ## get all the files match the pattern, that is one file only
                layer_to_add = arcpy.mapping.Layer(os.path.join(file_dir, file_to_add_path))
                arcpy.mapping.AddLayer(df, layer_to_add, "TOP")
                layer_to_update = arcpy.mapping.ListLayers(mxd, layer_to_add.name, df)[0]  ## redefine layer_to_update as an arcpy object
                arcpy.mapping.UpdateLayer(df, layer_to_update, source_lyr, True)  ## update symbology by using the one of source_lyr (has to be an arpy object and not just the path)

                ## Manually change symbology of classified raster
                if layer_to_update.symbologyType == "RASTER_CLASSIFIED":
                    layer_to_update.symbology.classBreakValues = params['class_break_values']
                    layer_to_update.symbology.classBreakLabels = params['class_break_labels']

                arcpy.RefreshTOC()
                arcpy.RefreshActiveView()

        elif file_dir == mask_folder:

            source_lyr_path = mask_lyr_path

            source_lyr = arcpy.mapping.Layer(source_lyr_path)
            search_key = "%s_%d.dat" % (mask_file_base_name, yr)  ## construct search key with variable parts replaced by *

            file_to_add_path = glob.glob(search_key)[0]  ## get all the files match the pattern, that is one file only
            layer_to_add = arcpy.mapping.Layer(os.path.join(file_dir, file_to_add_path))
            arcpy.mapping.AddLayer(df, layer_to_add, "TOP")
            layer_to_update = arcpy.mapping.ListLayers(mxd, layer_to_add.name, df)[0]  ## redefine layer_to_update as an arcpy object
            arcpy.mapping.UpdateLayer(df, layer_to_update, source_lyr, True)  ## update symbology by using the one of source_lyr (has to be an arpy object and not just the path)

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()
    

#%% SAVE MXD

print('Saving MXD...')

mxd.saveACopy(output_mxd)



#%% ZOOM ON BOOKMARKS AND EXPORT PNGs

if params['export_pngs']:

    print('Exporting PNGs...')
    
    ## Turn of all layers in list
    for layer in arcpy.mapping.ListLayers(mxd, '', df):
        layer.visible = False
    arcpy.RefreshActiveView()
    
    ## Loop over the bookmarks
    for bkmk in arcpy.mapping.ListBookmarks(mxd):  
        df.extent = bkmk.extent   ## set extent of dataframe to extent of bookmark

        ## Loop over year of interest
        for yr in params['years']:

            search_key = 'LC_Class_HMM*%d*' % (yr)
            layer = arcpy.mapping.ListLayers(mxd, search_key, df)[0]
            placement = "TOP"
            layer.visible = True

            search_key = '*%s*%d*' % (params['method'],  yr)
            ## Loop over the layers in df matching that year
            for layer in arcpy.mapping.ListLayers(mxd, search_key, df):
                placement = "BOTTOM"
                layer.visible = True
                arcpy.RefreshActiveView()
                file_name_out = "%s_%s.png" % (bkmk.name.split('_UTM')[0], os.path.splitext(layer.name)[0])   # build dynamic filename with the equivalent of sprintf
                arcpy.mapping.ExportToPNG(mxd, os.path.join(out_folder, file_name_out))  ## export mxd (zoomed on bookmark extent, and with only one layer available) to png
                layer.visible = False
                arcpy.RefreshActiveView()

print(toc(start_time))

