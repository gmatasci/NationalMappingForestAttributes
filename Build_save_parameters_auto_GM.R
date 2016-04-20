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
#  paramsGL$zones <- c('UTM9S', 'UTM12N', 'UTM13S', 'UTM20S')
 # paramsGL$zones <- c('UTM9S', 'UTM12N')
 # paramsGL$zones <- c('UTM13S')
 paramsGL$TARGET_YEAR <- 2010
 paramsGL$syr <- '2010'
 paramsGL$unit.size <- 25   ## use 75 to have Harold's original setting
 
 ##----------------------
 ## DIRECTORIES
 ##----------------------
 base_dir <- 'D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working'
 data_dir <- 'E:'
 LOP_dir <- file.path(data_dir,'LOP_data', fsep = .Platform$file.sep)
 Landsat_dir <- 'E:/NTEMS'
 topo_dir <- file.path(data_dir, 'TopoData', fsep = .Platform$file.sep)
 base_results_dir <- file.path(base_dir,'results',fsep = .Platform$file.sep)
 base_figures_dir <- file.path(base_dir,'figures',fsep = .Platform$file.sep)
 base_wkg_dir <- file.path(base_dir,'wkg', fsep = .Platform$file.sep)
 ##----------------------
 
 if (! file.exists(base_results_dir)){ dir.create(base_results_dir, showWarnings = F, recursive = T)}
 if (! file.exists(base_figures_dir)){ dir.create(base_figures_dir, showWarnings = F, recursive = T)}
 if (! file.exists(base_wkg_dir)){   dir.create(base_wkg_dir, showWarnings = F, recursive = T)}
 
 param_file = file.path(base_wkg_dir, 'AllUTMzones_paramsGL.Rdata', fsep = .Platform$file.sep) 
 
 save(paramsGL, base_dir, data_dir, LOP_dir, Landsat_dir, topo_dir, base_results_dir, base_figures_dir, base_wkg_dir, file = param_file)
 print(param_file)
 
 
 ## To run the entire workflow at once: not exactly working!
#  source('D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/prog0a_check_LiDARdata_row_duplicates.R')
#  source('D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/prog0_lidar_plots_sampling_auto_GM.R')
#  source('D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/prog1_training_validation_selection_auto_GM.R')
#  source('D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/prog2_lidar_plots_footprint_extract_exvars_auto_GM.R')
#  source('D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/prog3_rf_gnn_models_auto_GM.R')
#  source('D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/prog4_impute_the_whole_enclilada_auto_GM.R')
 
 