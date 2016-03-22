 ##########################################################
  # Project Name: CA_IMPUTATION #
  # Author: Geordie Hoabart    ghobart@nrcan.gc.ca         
  # File Name: build_save_aparams.R                             
  # Objective: A central repository for all global variables
  # Modifications                                          
  ##########################################################

 rm(list=ls())
 
 paramsGL <- list()
 paramsGL$zones <- c('UTM8N', 'UTM9N', 'UTM9S', 'UTM10N', 'UTM10S', 'UTM11N', 'UTM11S', 'UTM12N', 'UTM12S', 'UTM13S','UTM14S','UTM15S', 'UTM16S' ,'UTM17S', 'UTM18S' ,'UTM19S' ,'UTM20S' ,'UTM21S')
 # paramsGL$zones <- c('UTM12S', 'UTM13S', 'UTM14S')
 # paramsGL$zones <- c('UTM13S')
 # paramsGL$zones <- c('UTM9S', 'UTM12N')
 # paramsGL$zones <- c('UTM12N', 'UTM21S')
 paramsGL$TARGET_YEAR <- 2010
 paramsGL$syr <- '2010'
 paramsGL$unit.size <- 25   ## use 75 to have Harold's original setting
 
 ##----------------------
 ## DIRECTORIES
 ##----------------------
 base_dir = 'D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working'
 data_dir = file.path(base_dir,'data', fsep = .Platform$file.sep)
 LOP_dir = file.path(data_dir,'LOP_data', fsep = .Platform$file.sep)
 TC_dir = file.path(data_dir,'tc_composites', fsep = .Platform$file.sep)
 # change_dir = file.path(data_dir,'change', fsep = .Platform$file.sep)
 change_dir = 'E:/NTEMS'
 # topo_dir = file.path(data_dir, 'topo', fsep = .Platform$file.sep)  # Zald's  (JUST UTM13S)
 topo_dir = file.path(data_dir, 'TopoData', fsep = .Platform$file.sep)  # Geordie's  (ALL UTM ZONES)
 base_results_dir = file.path(base_dir,'results',fsep = .Platform$file.sep)
 base_wkg_dir = file.path(base_dir,'wkg', fsep = .Platform$file.sep)
 ##----------------------
 
 if (! file.exists(base_results_dir)){ dir.create(base_results_dir, showWarnings = F, recursive = T)}
 if (! file.exists(base_wkg_dir)){   dir.create(base_wkg_dir, showWarnings = F, recursive = T)}
 
 param_file = file.path(base_wkg_dir, 'AllUTMzones_paramsGL.Rdata', fsep = .Platform$file.sep) 
 
 save(paramsGL, base_dir, data_dir, LOP_dir, TC_dir, change_dir, topo_dir, base_results_dir, base_wkg_dir, file = param_file)
 print(param_file)