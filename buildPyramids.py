import arcpy, os, sys, multiprocessing
from arcpy import env

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*

tic = tic()

# Set the current workspace
#env.workspace = r'E:\NTEMS\UTM_results\UTM_13S_percentage_first_returns_above_2m'
#env.workspace = r'E:\NTEMS\UTM_results\RF_corrected'
env.workspace = r'E:\NTEMS\UTM_results\YAI_corrected'

## Build pyramids for all .dat files in folder
for layer in arcpy.ListFiles("*.dat"):
    
#    if layer.endswith("_stdev.dat"):   ## skip layers of standard deviation...
#        continue    ## ...by going to next round of loop
    
    print layer
    arcpy.BuildPyramids_management(os.path.join(env.workspace, layer), skip_existing="SKIP_EXISTING")  ## do not create pyramids if they already exist
	
toc(tic)
	
## Test with parallel approach (not exactly working!)
#pool = multiprocessing.Pool(24, None, None, 1)
#
#pool.map(arcpy.BuildPyramids_management, arcpy.ListFiles())
#
## closing pool and waiting for each task to be finished
#pool.close()
#pool.join()