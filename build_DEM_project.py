#########################################################
#  apply_543sybology2project_layers.py
#
#  The name says it all. To run just load your surface reflectance layers
#  into a arc project and save. Set the input mxd, output mxd and layer path below and run.
#
#  by Geordie Hobart
#     Forest Geomatics, Pacific Forestry Centre
#     Canadian Forest Service, Natural Resources Canada
#     286 - 506 Burnside Rd W.
#     Victoria, BC
#     V8Z 1M5
#     ghobart@nrcan.gc.ca
#     250-298-2403
#################################################################   


import arcpy.mapping, sys,os, glob, fnmatch

from arcpy import env
from arcpy.sa import *

TARGET_YEARS = range(1984,2013,1)

skip =[]

live = 0

finland = r'F:\Mikko'
rt_dirs= [r'F:\Mikko' ]



##input_mxd = os.path.join(finland, r'UTM_results\empty.mxd')
input_mxd = os.path.join(finland, r'layers\empty.mxd')
output_mxd =  os.path.join( rt_dirs[0], 'UTM35S_Proxy_change.mxd')
print input_mxd


nlyr = os.path.join(finland, r'layers\NewGroupLayer.lyr')
sym543_lyr = os.path.join(finland, r'layers\RGB_543_PC2_v3.lyr')
ch_mag_lyr = os.path.join(finland, r'layers\Change_magnitude_variation.dat.lyr')
ch_by_yr =os.path.join(finland, r'layers\Change_by_year.lyr')
ch_persis_lyr = os.path.join(finland, r'layers\Change_persistence.dat.lyr')
rate_lyr = os.path.join(finland, r'layers\Change_rate.dat.lyr')
grt_yr_lyr = os.path.join(finland, r'layers\Greatest_Change_Year.dat.lyr')
pst_ch_erate_lyr = os.path.join(finland, r'layers\PostChange_evolution_rate.dat.lyr')
pre_ch_erate_lyr = os.path.join(finland, r'layers\PreChange_evolution_rate.dat.lyr')
ch_year = os.path.join(finland, r'layers\Change_Year.dat.lyr')
sn_lyr = os.path.join(finland, r'layers\snow2.lyr')
a2p_lyr=os.path.join(finland, r'layers\area2process.lyr')
nd_lyr = os.path.join(finland, r'layers\No_data.lyr')

zones = {}
wkgMks = {}


i_srch_key = 'SRef_UTM*.dat'
dm_srch_key = 'dmask_UTM*.dat'
wkg_srch_key = 'utm*_wkg_data_mask.dat'
vld_srch_key = 'UTM*_vld_dmask.dat'
wtr_srch_key = 'UTM*_wtr_mask.dat'
s_srch_key = 'UTM*_snow_final.dat'
nd_srch_key = 'SRef_*_No_data_*_filter.dat'

s_srch_key2 = 'SRef_*_Snow_mask.dat'



mxd = arcpy.mapping.MapDocument(input_mxd)
df = arcpy.mapping.ListDataFrames(mxd)[0]
mxd.activeView = df.name
mxd.title = df.name

for rt_dir in rt_dirs:
        for cdir in os.listdir(rt_dir):
                tpath = os.path.join(rt_dir, cdir)
                if cdir.find('UTM') == 0 and os.path.isdir(tpath):
                        zn_num = cdir.split('_')[1]
                        
                        zn = cdir
                        
##                        if os.path.exists(output_mxd):
##                                continue
                       

                        print 'adding 2 ' +zn
                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = str(zn)
                        arcpy.mapping.AddLayer(df,  groupLayer, "BOTTOM")

                        chg_lyr ='ChangeResults'

                        ##print ' adding change'
                        ch_lyr = ''
                        tgt_grlyr = ''
                        UTMlyrHT = {}
                        for lyr in arcpy.mapping.ListLayers(mxd):
                                print "found" + lyr.name
                                yr = ''
                                if lyr.name.find(zn) == 0 :
                                       UTMlyrHT[lyr.name] = lyr
                                       tgt_grlyr = lyr
                                       print tgt_grlyr



                        chg_lyr_nm ='ChangeResults'
                        print ' adding change'
                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = chg_lyr_nm
                        arcpy.mapping.AddLayerToGroup(df, tgt_grlyr,  groupLayer, "TOP")

                        px_lyr_nm ='Proxy'
                        print ' adding proxy'
                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = px_lyr_nm
                        arcpy.mapping.AddLayerToGroup(df, tgt_grlyr,  groupLayer, "TOP")

                        sn_lyr_nm ='Snow'
                        print ' adding snow'
                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = sn_lyr_nm
                        arcpy.mapping.AddLayerToGroup(df, tgt_grlyr,  groupLayer, "TOP")

                        a2p = os.path.join(tpath,'area_to_process','UTM'+zn_num+'_area_to_process_0123.dat')
                        arcpy.BuildPyramidsandStatistics_management(a2p,"INCLUDE_SUBDIRECTORIES","BUILD_PYRAMIDS",\
                                                                                    "CALCULATE_STATISTICS","NONE","#","NONE","1","1","#",\
                                                                                    "-1","NONE","NEAREST","NONE","75","SKIP_EXISTING")                           
                        arcpy.BuildRasterAttributeTable_management(a2p) 
                        result = arcpy.MakeRasterLayer_management(a2p, 'Area_to_process')
                        layer = result.getOutput(0)
                        arcpy.ApplySymbologyFromLayer_management(layer, a2p_lyr )
                        arcpy.mapping.AddLayer(df , layer, 'TOP')
                        arcpy.mapping.AddLayerToGroup(df, tgt_grlyr,  layer, "TOP")

                        tgt_grlyr = ''
                        sn_grlyr=''
                        for lyr in arcpy.mapping.ListLayers(mxd):
                                if lyr.name == chg_lyr:
                                        ch_grlyr = lyr
                                if lyr.name == px_lyr_nm:
                                        px_grlyr = lyr       
                                if lyr.name == sn_lyr_nm:
                                        sn_grlyr = lyr 

                        ##for zn in UTMlyrHT.keys():
                        ##        tgt_grlyr = UTMlyrHT[zn]
                        for yr in TARGET_YEARS:
                                print yr
                                groupLayer = arcpy.mapping.Layer(nlyr)
                                groupLayer.name = str(yr)
                                arcpy.mapping.AddLayerToGroup(df, px_grlyr,  groupLayer, "BOTTOM")

                              
                                                        
                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = 'PreChange'
                        arcpy.mapping.AddLayerToGroup(df, ch_grlyr,  groupLayer, "TOP")

                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = 'PostChange'
                        arcpy.mapping.AddLayerToGroup(df, ch_grlyr,  groupLayer, "BOTTOM")

                        groupLayer = arcpy.mapping.Layer(nlyr)
                        groupLayer.name = 'Change'
                        arcpy.mapping.AddLayerToGroup(df, ch_grlyr,  groupLayer, "TOP")

                        chHT = {}
                        YRlyrHT = {}
                        for lyr in arcpy.mapping.ListLayers(mxd):
                                print lyr.name
                                if lyr.name.find('UTM') == 0:
                                        znName = lyr.name
                                        YRlyrHT[znName] = {}
                                if lyr.name.find('19') == 0 or lyr.name.find('20') == 0:
                                        YRlyrHT[znName][lyr.name] = lyr
                                if lyr.name.find('Change') > -1:
                                        chHT[lyr.name] = lyr        




                        cmbnd_change = os.path.join(tpath,'Results','SRef_'+zn_num+'_Combined_changes_Multiyear.dat')
                        if not os.path.exists(cmbnd_change):
                                        cmbnd_change = os.path.join(tpath,'Results','Combined_changes_Multiyear.dat')
                        print cmbnd_change
                                                   

                        arcpy.BuildPyramidsandStatistics_management(cmbnd_change,"INCLUDE_SUBDIRECTORIES","BUILD_PYRAMIDS",\
                                                                                    "CALCULATE_STATISTICS","NONE","#","NONE","1","1","#",\
                                                                                    "-1","NONE","NEAREST","NONE","75","SKIP_EXISTING")

                                       

                        tgt_dir = os.path.join(tpath,'Results','proxy_values')
                        os.chdir(tgt_dir)
                        res = glob.glob(i_srch_key)

                        nd_dir = os.path.join(tpath,'no_data')
                        os.chdir(nd_dir)
                        nd_res = glob.glob(nd_srch_key)
                        ndHt = {}
                        for nd in nd_res:
                                #SRef_35S_No_data_1984_filter
                                yr = nd.split('_')[4]
                                ndHt[yr] = os.path.join(nd_dir, nd)
                                
                 
                        
                        for rs in res:
                                        print rs
                                        yr = rs.split('_')[2]
                        ##                print
                                        
                                        grpLyr = YRlyrHT[zn][yr]
                                        print grpLyr
                                        data = os.path.join(tgt_dir,rs)
                                        arcpy.BuildPyramidsandStatistics_management(data,"INCLUDE_SUBDIRECTORIES","BUILD_PYRAMIDS",\
                                                                                 "CALCULATE_STATISTICS","NONE","#","NONE","1","1","#",\
                                                                                    "-1","NONE","NEAREST","NONE","75","SKIP_EXISTING")
                                
                                        print data
                                        result = arcpy.MakeRasterLayer_management(data, rs[:-4])
                                        layer = result.getOutput(0)
                                        arcpy.ApplySymbologyFromLayer_management(layer, sym543_lyr )
                                        arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'BOTTOM')

                                        # add chage layer

                                        
                                        cmbnd_change_band = os.path.join(cmbnd_change, yr)
                        ##                print cmbnd_change_band
                                        result = arcpy.MakeRasterLayer_management(cmbnd_change_band, 'Change_'+yr)
                                        layer = result.getOutput(0)
                                        layer.visible=False
                        ##                print "'"+ch_mag_lyr+"'"
                                        arcpy.ApplySymbologyFromLayer_management(layer, ch_by_yr )
                                        arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'TOP')                


                                        nd = ndHt[yr]
                                        result = arcpy.MakeRasterLayer_management(nd, 'no_data')
                                        layer = result.getOutput(0)
                                        arcpy.ApplySymbologyFromLayer_management(layer, nd_lyr )
                                        arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'TOP')





                        tgt_dir = os.path.join(tpath,'Snow')
                        os.chdir(tgt_dir)
                        res = glob.glob(s_srch_key2)

                        prev_yr_snow =''

                        grpLyr = sn_grlyr
                        data = os.path.join(tgt_dir,res[0])
                        arcpy.BuildPyramidsandStatistics_management(data,"INCLUDE_SUBDIRECTORIES","BUILD_PYRAMIDS",\
                                                                                    "CALCULATE_STATISTICS","NONE","#","NONE","1","1","#",\
                                                                                    "-1","NONE","NEAREST","NONE","75","SKIP_EXISTING")
                        print data
                        result = arcpy.MakeRasterLayer_management(data, rs[:-4])
                        layer = result.getOutput(0)
                        arcpy.ApplySymbologyFromLayer_management(layer, sn_lyr )
                        arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'TOP')

                        ##sn_grlyr
                                               
                        ##for rs in res:
                        ##                print rs
                        ##                yr = rs.split('_')[1]
                        ####                UTM10S_1984_snow_final
                        ##                
                        ##                grpLyr = YRlyrHT[zn][yr]
                        ##                print grpLyr
                        ##                data = os.path.join(tgt_dir,rs)
                        ##                arcpy.BuildPyramidsandStatistics_management(data,"INCLUDE_SUBDIRECTORIES","BUILD_PYRAMIDS",\
                        ##                                                            "CALCULATE_STATISTICS","NONE","#","NONE","1","1","#",\
                        ##                                                            "-1","NONE","NEAREST","NONE","75","SKIP_EXISTING")
                        ##        
                        ##                print data
                        ##                result = arcpy.MakeRasterLayer_management(data, rs[:-4])
                        ##                layer = result.getOutput(0)
                        ##                arcpy.ApplySymbologyFromLayer_management(layer, sn_lyr )
                        ##                arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'TOP')
                        ##
                        ##                # add chage layer
                        ##
                        ##                if prev_yr_snow <> '' :
                        ##                        
                        ##                        result = arcpy.MakeRasterLayer_management(prev_yr_snow, prev_name)
                        ##                        layer = result.getOutput(0)
                        ##        ##              
                        ##                        arcpy.ApplySymbologyFromLayer_management(layer, sn_lyr )
                        ##                        arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'TOP')
                        ##                prev_yr_snow = data
                        ##                prev_name = rs[:-4]





                        tgt_dir = os.path.join(tpath,'Results','Change_metrics')
                        os.chdir(tgt_dir)
                        res = glob.glob('*.dat')

                        for rs in res:
                                data = os.path.join(tgt_dir,rs)
                                # define the group layer
                                grpLyr =  chHT['Change']
                                
                                if rs.find('PreChange')> -1:
                                       grpLyr =  chHT['PreChange']
                                if rs.find('PostChange')> -1:
                                       grpLyr =  chHT['PostChange']
                                
                                        
                                 # define the symbology layer       
                                sym_lyr = ch_persis_lyr
                                if rs.find('magnitude_variation')> -1:
                                       sym_lyr =  ch_mag_lyr
                                if rs.find('rate')> -1:
                                        sym_lyr =  rate_lyr
                                        if rs.find('PostChange_evolution')> -1:
                                               sym_lyr =  pst_ch_erate_lyr
                                        if rs.find('PreChange_evolution')> -1:
                                               sym_lyr =  pre_ch_erate_lyr
                                if rs.find('Year')> -1:
                                        arcpy.BuildRasterAttributeTable_management(data) 
                                        sym_lyr =  grt_yr_lyr
                        ##        if rs.find('Change_persistence')> -1:
                        ##                arcpy.BuildRasterAttributeTable_management(data) 
                        ##                sym_lyr =  ch_persis_lyr                # add the layer       
                                data = os.path.join(tgt_dir,rs)
                                arcpy.BuildPyramidsandStatistics_management(data,"INCLUDE_SUBDIRECTORIES","BUILD_PYRAMIDS",\
                                                                                "CALCULATE_STATISTICS","NONE","#","NONE","1","1","#",\
                                                                                "-1","NONE","NEAREST","NONE","75","SKIP_EXISTING")
                                print data
                                result = arcpy.MakeRasterLayer_management(data, rs[:-4])
                                layer = result.getOutput(0)
                                arcpy.ApplySymbologyFromLayer_management(layer, sym_lyr )
                                arcpy.mapping.AddLayerToGroup(df , grpLyr, layer, 'BOTTOM')

                        print 'saving'
mxd.saveACopy(output_mxd)

