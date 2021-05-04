# NationalMappingForestAttributes
R and Python scripts for the Postdoc project (collaboration between UBC's Integrated Remote Sensing Studio and the Canadian Forest Service).

## Goal
Develop machine learning models (Random Forests, imputation methods) to relate time-series of Landsat images and LiDAR data to ultimately map forest attributes across Canada over a period of 30+ years, improving forest monitoring and management practices.

## Approach and main results
### Single year map
<img src="https://github.com/gmatasci/NationalMappingForestAttributes/blob/master/Figures/NatMappingForeestAttr_Workflow_wLiDARtransectFieldPlots.png" width="800">

<img src="https://github.com/gmatasci/NationalMappingForestAttributes/blob/master/Figures/BorealMaps.png" width="800">

### Multitemporal extension
<img src="https://github.com/gmatasci/NationalMappingForestAttributes/blob/master/Figures/TemporalMapping_allData_trendFitting.png" width="800">

<img src="https://github.com/gmatasci/NationalMappingForestAttributes/blob/master/Figures/Fig7.png" width="800">

## Script description
- build_BorealMaps_YAI_2010.py: Build MXD with the 10 forest attributes maps for 2010 masked on the boreal

- build_Pyramids_UTM_maps.py: Build pyramids for all forest attributes maps in each UTM zone saved in 'E:\NTEMS\UTM_results'

- Build_save_parameters.R: Set global parameters and base directories for the suite of R scripts starting with 'prog0_lidar_plots_sampling_auto_GM.R'

- exportMaps_YAI_2010.py: Export to PNG the 2010 maps of some specified attributes of interest zoomed to a given bookmark 

- exportTimeSeries_CAN_YAI_1984_2012.py: Export to PNG the time-series of specified attributes of interest zoomed on bookmarks across Canada (maps of 2nd paper)

- exportTimeSeries_elevp95_RF_1984_2012_YAI_2010.py: Export to PNG the time-series of some specified attributes of interest zoomed to a given bookmark 

- Functions_NatMapping_Python.py: Define Python functions to make available to all py scripts

- Functions_NatMapping_R.R: Define R functions to make available to all R scripts

- merge_transects.py: Merge into a single shp all the individual shps containing the transects in each UTM zone

- prog0_lidar_plots_sampling.R: Define sampling grid over LiDAR transect

- prog0a_check_LiDARdata_row_duplicates.R: Check for row duplicates in all LiDAR csv files for all UTM zones

- prog0b_create_shp_from_csv.R: Create shp from CSV files and save new CSV files (lidar and forest attributes) with new unique_id for new LiDAR data from BC

- prog1_training_validation_selection.R: LiDAR training and validation plot selection

- prog1a_volume_biomass_check_recompute.R: check the issue with incoherence of Vol and Biomass values between Boreal/Non-boreal ecozones. Wrong equations have been used and this script recomputes the correct values on the transect for these 2 attributes.


- prog2_lidar_plots_footprint_extract_exvars.R: Extract values for all explanatory rasters spatially coincident with sampled lidar plot polygons for both training and validation plots


- prog3_model_selection.R: Descriptive statistics and plots, variable selection and Random Forest prediction/imputation

- prog4a_map_donor_plot_distance.py: Computes the distance from each pixel to its donor plot (in Lat/Long degree units) and also saves Lat and Long of donor plot

- prog4_check_maps.R: Check that predicted values on maps match those predicted on plots

- prog5a_rename_maps.py: renames maps after correction of biomass and volume values

- prog5_ecozone_stats.py: Compute descriptive stats by ecozone over the entire boreal for all the mapped attributes

- prog5_ecozone_stats_plots.R: Produce boxplots based on the descriptive stats extracted with 'prog5_ecozone_stats.py'

- prog6a_predict_on_random_pix.R: Predict forest attributes through time for randomly sampled pixels whose predictors are provided in a CSV file.

- prog6b_predict_on_val_pix.R: Predict forest attributes through time for validation pixels whose predictors are provided in a CSV file.

- prog6_temporal_analysis.R: Plot time-series of mapped attributes with boxplots of an ensemble of pixels sampled over large areas or with lines for individual pixel
