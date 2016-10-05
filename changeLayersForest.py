# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 17:17:49 2016

@author: gmatasci
"""

#import arcpy
#from arcpy import env
#
## Set the current workspace
#env.workspace = "path/to/workspace"
#
## Set layer to apply symbology to
#inputLayers = ["in_layer_first.lyr","in_layer_second.lyr","in_layer_third.lyr"]
#
## Set layer that output symbology will be based on
#symbologyLayer = "utm10n_fsum.lyr"
#
## Apply the symbology from the symbology layer to the input layer
#for layer in inputLayers:
#    arcpy.ApplySymbologyFromLayer_management (layer, symbologyLayer)

import arcpy

#Get the current Map Document
mxd = arcpy.mapping.MapDocument("CURRENT")

# Script arguments
Template_Layer = arcpy.GetParameterAsText(0)
LayerList = arcpy.GetParameterAsText(1)
Layers_to_Symbolize = LayerList.split(";")

# Process: Apply Symbology From Layer
for UpdateLayer in Layers_to_Symbolize:
    arcpy.AddMessage("Updating: " + UpdateLayer)
    arcpy.ApplySymbologyFromLayer_management(UpdateLayer,Template_Layer)

# Refresh the Table of Contents to reflect the change
arcpy.RefreshTOC()

#Delete the MXD from memory
del mxd