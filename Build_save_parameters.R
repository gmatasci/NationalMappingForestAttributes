## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hobart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: Build_save_parameters.R                         
## Objective: A central repository for all global variables

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - check final UTM zones to sample (now 17 with the removal of 11S)

## SOLVED:

#### GLOBAL PARAMETERS -------------------------------------------------------

print('Build and save parameters')

rm(list=ls())

paramsGL <- list()
paramsGL$zones <- c('UTM8N', 'UTM9N', 'UTM9S', 'UTM10N', 'UTM10S', 'UTM11N', 'UTM12N', 'UTM12S', 'UTM13S','UTM14S','UTM15S', 'UTM16S' ,'UTM17S', 'UTM18S' ,'UTM19S' ,'UTM20S' ,'UTM21S')
# paramsGL$zones <- c('UTM10N', 'UTM10S')
paramsGL$TARGET_YEAR <- 2010
paramsGL$global.seed <- 2010  ## to have the same initialization for any random process
paramsGL$unit.size <- 25   ## use 75 to have Harold's original setting

#### DIRECTORIES -------------------------------------------------------------

base_dir <- 'D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR'
data_dir <- 'E:'
LOP_dir <- file.path(data_dir,'LOP_data', fsep = .Platform$file.sep)
Landsat_dir <- 'E:/NTEMS'
topo_dir <- file.path(data_dir, 'TopoData_v2', fsep = .Platform$file.sep)
base_results_dir <- file.path(base_dir,'results',fsep = .Platform$file.sep)
base_figures_dir <- file.path(base_dir,'figures',fsep = .Platform$file.sep)
base_wkg_dir <- file.path(base_dir,'wkg', fsep = .Platform$file.sep)

if (! file.exists(base_results_dir)){ dir.create(base_results_dir, showWarnings = F, recursive = T) }
if (! file.exists(base_figures_dir)){ dir.create(base_figures_dir, showWarnings = F, recursive = T) }
if (! file.exists(base_wkg_dir)){ dir.create(base_wkg_dir, showWarnings = F, recursive = T) }

#### SAVE PARAMETERS ---------------------------------------------------------

param_file = file.path(base_wkg_dir, 'AllUTMzones_paramsGL.Rdata', fsep = .Platform$file.sep) 
save(paramsGL, base_dir, data_dir, LOP_dir, Landsat_dir, topo_dir, base_results_dir, base_figures_dir, base_wkg_dir, file = param_file)
print(param_file)

