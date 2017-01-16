"""
Project Name: NationalMappingForestAttributes
Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
File Name: exportTimeSeries_elevp95_RF_1984_2012_YAI_2010.py                          
Objective: Export to PNG the time-series of some specified attributes of interest zoomed to a given bookmark 
"""

## IMPORT MODULES ---------------------------------------------------------------

import arcpy, os, sys, glob
import numpy as np

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*

## PARAMETERS  ------------------------------------------------

mxds_folder = r"D:\Research\ANALYSES\NationalMappingForestAttributes\WKG_DIR_NationalMappingForestAttributes\MXDs"

root_dir = os.path.join(mxds_folder, "TimeSeries", "UTM13S_Elev_p95_timeseries_1984_2012")

input_mxd = os.path.join(root_dir, r'empty.mxd')  ## path to empty mxd that has to have a scale bar (will adapt automatically) and bookmarks already set
output_mxd =  os.path.join(root_dir, r"UTM13S_Elev_p95_timeseries_1984_2012.mxd")  ## path to final mxd to save results

params = {}

params['years'] = np.arange(1984, 2013, 1)
#params['years'] = np.arange(1984, 1986, 1)

params['symbol_type'] = "continuous" 
#params['symbol_type'] = "classified" 

params['class_break_values'] = [2000, 7000, 10000, 13000, 16000, 60000]
params['class_break_labels'] = ["2 to 7 m", "7 to 10 m", "10 to 13 m", "13 to 16 m", "> 16 m"]

#params['targ'] = "cover_2m"
params['targ'] = "elev_p95"

params['methods'] = ["RF", "YAI"]

#params['export_pngs'] = False
params['export_pngs'] = True

in_folder = r"E:\NTEMS\UTM_results\UTM_13S"

in_folder_RF = os.path.join(in_folder, r'RF') 
in_file_base_name_RF = "UTM_13S_"+params['targ']+"_RF"
in_folder_YAI = os.path.join(in_folder, r'YAI') 
in_file_base_name_YAI = "UTM_13S_"+params['targ']+"_YAI"

mask_folder = r"E:\NTEMS\UTM_results\VLCE\13S_HMM"
mask_file_base_name = "LC_Class_HMM_13S"

lyrs_folder = os.path.join(mxds_folder, "LYRs")

elev_p95_lyr_contin_path = os.path.join(lyrs_folder, r'Elev_p95_contin_5to19m.lyr')  ## path to layer file providing the symbology for the RGB composite
elev_p95_lyr_class_path = os.path.join(lyrs_folder, r'Elev_p95_class_0to23m.lyr')  ## path to layer file providing the symbology for the RGB composite
mask_lyr = os.path.join(lyrs_folder, r'Mask_4_TreedClasses.lyr')  ## path to layer file providing the symbology for the categorical rasters (unique values)

out_folder = os.path.join(root_dir, r"Test_nr1")


## START --------------------------------------------------------------------------

print('Start')

tic = tic()

file_dirs = [in_folder_RF, in_folder_YAI, mask_folder]
in_file_base_names = [in_file_base_name_RF, in_file_base_name_YAI, mask_file_base_name]

## Create directory for output .png files
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

## Get MXD and dataframe
mxd = arcpy.mapping.MapDocument(input_mxd)
df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]  ## assuming there is only 1 df (called Layer, the default)
mxd.activeView = df.name
mxd.title = df.name

## LOAD LAYERS AND APPLY SYMBOLOGY ----------------------------------------------

print('Loading layers and applying symbology...')

if params['symbol_type'] == "continuous":
    elev_p95_lyr_path = elev_p95_lyr_contin_path
elif params['symbol_type'] == "classified": 
    elev_p95_lyr_path = elev_p95_lyr_class_path

for idx_dir, file_dir in enumerate(file_dirs):  ## iterate over the directories
    
    arcpy.env.workspace = file_dir   ## set this directory as workspace on which arcpy.ListFiles() will be called on
    os.chdir(file_dir)
    if file_dir == in_folder_RF:  
        years_to_add = params['years']
        source_lyr_path = elev_p95_lyr_path     
    elif file_dir == in_folder_YAI:
        years_to_add = [2010]
        source_lyr_path = elev_p95_lyr_path    
    elif file_dir == mask_folder:
        years_to_add = params['years']
        source_lyr_path = mask_lyr     
    
    source_lyr = arcpy.mapping.Layer(source_lyr_path)
 
    for yr in years_to_add:  
        
        search_key = in_file_base_names[idx_dir]+'_%d.dat' % yr  ## construct search key with variable parts replaced by *
        file_to_add_path = glob.glob(search_key)[0]  ## get all the files match the pattern, that is one file only
        layer_to_add = arcpy.mapping.Layer(os.path.join(file_dir, file_to_add_path))
        arcpy.mapping.AddLayer(df, layer_to_add, "TOP")
        layer_to_update = arcpy.mapping.ListLayers(mxd, layer_to_add.name, df)[0]   ## redefine layer_to_update as an arcpy object 
        arcpy.mapping.UpdateLayer(df, layer_to_update, source_lyr, True)  ## update symbology by using the one of source_lyr (has to be an arpy object and not just the path)
        
        ## Manually change symbology of classified raster
        if layer_to_update.symbologyType == "RASTER_CLASSIFIED":
            layer_to_update.symbology.classBreakValues = params['class_break_values']
            layer_to_update.symbology.classBreakLabels = params['class_break_labels']
        
        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
    

##  SAVE MXD -------------------------------------------------------------

print('Saving MXD...')

mxd.saveACopy(output_mxd)

## ZOOM ON BOOKMARKS AND EXPORT PNGs --------------------------------------

if params['export_pngs']:

    print('Exporting PNGs...')
    
    ## Turn of all layers in list
    for layer in arcpy.mapping.ListLayers(mxd, '', df):
        layer.visible = False
    arcpy.RefreshActiveView()
    
    ## Loop over the bookmarks
    for bkmk in arcpy.mapping.ListBookmarks(mxd):  
        df.extent = bkmk.extent   ## set extent of dataframe to extent of bookmark
        
        for method in params['methods']:
        
            ## Loop over year of interest
            for yr in params['years']:     
                
                ## Only valid layers are all year for RF or 2010 for YAI
                if method == 'RF' or (method == 'YAI' and yr == 2010):
                    
                    search_key = '*%d*' % (yr)  
                    
                    ## Loop over the layers in df matching that year
                    for layer in arcpy.mapping.ListLayers(mxd, search_key, df): 
                        if "LC_Class_HMM" in layer.name:
                            placement = "TOP"
                            layer.visible = True   
                        elif method in layer.name:
                            placement = "BOTTOM"
                            layer.visible = True   
                        else:
                            layer.visible = False
                        arcpy.RefreshActiveView()   
                    file_name_out = "%s_%s_%d.png" % (bkmk.name, method, yr)   # build dynamic filename with the equivalent of sprintf 
                    arcpy.mapping.ExportToPNG(mxd, os.path.join(out_folder, file_name_out))  ## export mxd (zoomed on bookmark extent, and with only one layer available) to png
            
                    for layer in arcpy.mapping.ListLayers(mxd, '', df):
                        layer.visible = False
                        arcpy.RefreshActiveView()

toc(tic)


