- build_BorealMaps_YAI_2010.py: Build MXD with the 10 forest attributes maps for 2010 masked on the boreal

- build_Pyramids_UTM_maps.py: Build pyramids for all forest attributes maps in each UTM zone saved in 'E:\NTEMS\UTM_results'

- Build_save_parameters.R: Set global parameters and base directories for the suite of R scripts starting with 'prog0_lidar_plots_sampling_auto_GM.R'

- exportMaps_YAI_2010.py: Export to PNG the 2010 maps of some specified attributes of interest zoomed to a given bookmark 

- exportTimeSeries_elevp95_RF_1984_2012_YAI_2010.py: Export to PNG the time-series of some specified attributes of interest zoomed to a given bookmark 

- Functions_NatMapping_Python.py: Define Python functions to make available to all py scripts

- Functions_NatMapping_R.R: Define R functions to make available to all R scripts

- mergeTransects.py: Merge into a single shp all the individual shps containing the transects in each UTM zone

- prog0_lidar_plots_sampling.R: Define sampling grid over LiDAR transect

- prog0a_check_LiDARdata_row_duplicates.R: Check for row duplicates in all LiDAR csv files for all UTM zones

- prog1_training_validation_selection.R: LiDAR training and validation plot selection

- prog2_lidar_plots_footprint_extract_exvars.R: Extract values for all explanatory rasters spatially coincident with sampled lidar plot polygons for both training and validation plots

- prog3_model_selection.R: Descriptive statistics and plots, variable selection and Random Forest prediction/imputation

- prog4_check_maps.R: Check that predicted values on maps match those predicted on plots

- prog5_ecozone_stats.py: Compute descriptive stats by ecozone over the entire boreal for all the mapped attributes

- prog5_ecozone_stats_plots.R: Produce boxplots based on the descriptive stats extracted with 'prog5_ecozone_stats.py'

- prog6_temporal_analysis.R: Plot time-series of mapped attributes with boxplots of an ensemble of pixels sampled over large areas or with lines for individual pixel

- prog6a_predict_on_random_pix.R: Predict forest attributes through time for randomly sampled pixels whose predictors are provided in a CSV file.

- prog6b_predict_on_val_pix.R: Predict forest attributes through time for validation pixels whose predictors are provided in a CSV file.