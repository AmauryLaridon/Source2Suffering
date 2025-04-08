#!/usr/bin/env python3                                                                                       #
# -----------------------------------------------------------------------------------------------------------#
# Main script to postprocess and visualise lifetime exposure to climate extremes data                        #
#                                                                                                            #
# Python translation of the MATLAB scripts of Thiery et al. (2021)                                           #
# The translation of the scripts has been first performed by Luke Grant for Grant et al.(2025)               #
# The modification of Luke Grant's scripts have been performed by Amaury Laridon for Laridon et al.(2025)    #
# -----------------------------------------------------------------------------------------------------------#
# Associated papers                                                                                          #
# -----------------------------------------------------------------------------------------------------------#
#                                                                                                            #
# - Thiery et al.(2021) - https://www.science.org/doi/abs/10.1126/science.abi7339                            #
# - Grant et al.(2025) - in review                                                                           #
# - Laridon et al.(2025) - in prep                                                                           #
#                                                                                                            #
# -----------------------------------------------------------------------------------------------------------#
# Summary and notes                                                                                          #
# -----------------------------------------------------------------------------------------------------------#
#                                                                                                            #
# Data types are defined in the variable names starting with:                                                #
#     df_     : DataFrame    (pandas)                                                                        #
#     gdf_    : GeoDataFrame (geopandas)                                                                     #
#     da_     : DataArray    (xarray)                                                                        #
#     d_      : dictionary                                                                                   #
#     sf_     : shapefile                                                                                    #
#     ...dir  : directory                                                                                    #
# -----------------------------------------------------------------------------------------------------------#
#%%----------------------------------------------------------------------------------------------------------#
# Libraries                                                                                                  #
#------------------------------------------------------------------------------------------------------------#

import xarray as xr
import pickle as pk
import time
import matplotlib as mpl
import mapclassify as mc
from copy import deepcopy as cp
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy as cr
import geopandas as gpd
import seaborn as sns 
import datetime
import time
import sys
import subprocess
import shutil

#%%----------------------------------------------------------------------------------------------------------#
# Flags - Define the Configuration of the Framework                                                          #
#------------------------------------------------------------------------------------------------------------#

# Global flags use to configure the framework

global flags

flags = {}

# Boolean to produce the exact configuration of the framework based on the different papers and methodologies 

global Thiery_2021, Grant_2025, Laridon_2025, Source2Suffering

Thiery_2021 = False 
Grant_2025 = False      
Laridon_2025 = False
Source2Suffering = True
env_value_paper= 0 

#-------------------------------------------------------------------------------------#
# Submit jobs that execute main.py with specific configurations - use on HPC only     #
#-------------------------------------------------------------------------------------#

# Use for specific paper's configuration jobs
env_value_paper = os.getenv("CONFIG_PAPER_VALUE")
if env_value_paper:
    Thiery_2021 = False 
    Grant_2025 = True 
    Laridon_2025 = False
    Source2Suffering = True 
    if env_value_paper == "thiery_2021":
        Thiery_2021 = True
    if env_value_paper == "grant_2025":
        Grant_2025 = True
    if env_value_paper == "laridon_2025":
        Laridon_2025 = True 
    if env_value_paper == "source2suffering":
        Source2Suffering = True   

#%%------------------------------------------------------------------------------------#
# Configuration of the Framework to reproduce the papers                               #
#--------------------------------------------------------------------------------------#
if env_value_paper:

    if Thiery_2021==True:
        flags['extr'] = 'heatwavedarea'
        flags['gmt'] = 'original'
        flags['rm'] = 'no_rm'
        flags['version'] = 'pickles'
        flags['run'] = 1
        flags['mask'] = 1
        flags['lifetime_exposure'] = 1
        flags['lifetime_exposure_pic'] = 1
        flags['emergence'] = 0
        flags['birthyear_emergence'] = 0
        flags['gridscale'] = 0
        flags['gridscale_le_test'] = 0
        flags['gridscale_country_subset'] = 0
        flags['global_emergence_recollect'] = 0
        flags['global_avg_emergence'] = 0
        flags['gdp_deprivation'] = 0
        flags['vulnerability'] = 0
        flags['plots'] = 0
        flags['reporting'] = 0

    if Grant_2025==True:
        flags['extr'] = 'all'
        flags['gmt'] = 'ar6_new'
        flags['rm'] = 'rm'
        flags['version'] = 'pickles_v3'
        flags['run'] = 0
        flags['mask'] = 0
        flags['lifetime_exposure'] = 0
        flags['lifetime_exposure_pic'] = 0
        flags['emergence'] = 1
        flags['birthyear_emergence'] = 0
        flags['gridscale'] = 1
        flags['gridscale_le_test'] = 0
        flags['gridscale_country_subset'] = 0
        flags['global_emergence_recollect'] = 0
        flags['global_avg_emergence'] = 0
        flags['gdp_deprivation'] = 0
        flags['vulnerability'] = 0
        flags['plots'] = 0
        flags['reporting'] = 0

    if Laridon_2025==True:
        print('Configuration of the Framework for Laridon et al.(2025) not settle going to Manual Configuration')
    
    if Source2Suffering==True:
        print('Configuration of the Framework for Laridon et al.(2025) not settle going to Manual Configuration')

#%%------------------------------------------------------------------------------------#
# Flags - Manual Configuration of the Framework                                        #
#--------------------------------------------------------------------------------------#

if not env_value_paper:

    #--------------------- Configuration ----------------------#

    flags['extr'] = 'heatwavedarea'                 # 0: all
                                                    # 1: burntarea
                                                    # 2: cropfailedarea
                                                    # 3: driedarea
                                                    # 4: floodedarea
                                                    # 5: heatwavedarea
                                                    # 6: tropicalcyclonedarea

    flags['gmt'] = 'ar6_new'                        # original: use Wim's stylized trajectory approach with max trajectory a linear increase to 3.5 deg                               
                                                    # ar6: substitute the linear max wth the highest IASA c7 scenario (increasing to ~4.0), new lower bound, and new 1.5, 2.0, NDC (2.8), 3.0
                                                    # ar6_new: works off ar6, but ensures only 1.5-3.5 with perfect intervals of 0.1 degrees (less proc time and data volume)

    flags['rm'] = 'rm'                              # no_rm: no smoothing of RCP GMTs before mapping
                                                    # rm: 21-year rolling mean on RCP GMTs
    
    flags['version'] = 'pickles_v3'                 # pickles: original version, submitted to Nature
                                                        # inconsistent GMT steps (not perfect 0.1 degree intervals)
                                                        # GMT steps ranging 1-4 (although study only shows ~1.5-3.5, so runs are inefficient)
                                                        # only 99.99% percentile for PIC threshold
                                                    # pickles_v2: version generated after submission to Nature in preparation for criticism/review
                                                        # steps fixed in load_manip to be only 1.5-3.5, with clean 0.1 degree intervals
                                                        # 5 percentiles for PIC threshold and emergence for each
                                                    # pickles_v3: version generated after the 2021 toolchains were taken away from hydra. could not longer use old pickles effectively

    flags['run'] = 0                                # 0: do not process ISIMIP runs (i.e. load runs pickle)
                                                    # 1: process ISIMIP runs (i.e. produce and save runs as pickle)

    flags['mask'] = 0                               # 0: do not process country data (i.e. load masks pickle)
                                                    # 1: process country data (i.e. produce and save masks as pickle)

    #--------- Thiery et al.(2021) Lifetime Exposure ----------#
    
    flags['lifetime_exposure'] = 0                  # 0: do not process ISIMIP runs to compute exposure across cohorts (i.e. load exposure pickle)
                                                    # 1: process ISIMIP runs to compute exposure across cohorts (i.e. produce and save exposure as pickle)   
                                                                        
    flags['lifetime_exposure_pic'] = 0              # 0: do not process ISIMIP runs to compute picontrol exposure (i.e. load exposure pickle)
                                                    # 1: process ISIMIP runs to compute picontrol exposure (i.e. produce and save exposure as pickle)

    #--------- Grant et al.(2025) Emergence of ULE ------------#
    
    flags['emergence'] = 0                          # 0: do not process ISIMIP runs to compute cohort emergence of unprecedentend lifetime exposure (i.e. load cohort exposure pickle)
                                                    # 1: process ISIMIP runs to compute cohort emergence of unprecedentend lifetime exposure (i.e. produce and save exposure as pickle)

    #--------------------- Produce Error ----------------------#

    flags['birthyear_emergence'] = 0                # 0: only run calc_birthyear_align with birth years from 1960-2020
                                                    # 1: run calc_birthyear_align with birth years from 1960-2100. Produces an error.  
                                                
    #--------- Grant et al.(2025) Gridscale analysis ----------#

                            
    flags['gridscale'] = 0                          # 0: do not process 1d scale analysis, load pickles
                                                    # 1: process grid scale analysis

    flags['gridscale_le_test'] = 0                  # 0: do not process the grid scale analysis testing diff versions of constant life expectancy
                                                    # 1: process grid scale analysis testing diff versions of constant life expectancy    
                                                                
    flags['gridscale_country_subset'] = 0           # 0: run gridscale analysis on all countries
                                                    # 1: run gridscale analysis on subset of countries determined in "get_gridscale_regions" and settings countries. Can only work if flags['gridscale'] = 1

    #------------- Grant et al.(2025) analysis ----------------#
    #--------------------- Produce Error ----------------------#

    flags['global_emergence_recollect'] = 0         # 0: do not load pickles of global emergence masks used for vulnerability assessment
                                                    # 1: load pickles                  

    flags['global_avg_emergence'] = 0               # 0: do not run d_global_emergence used for SI Figures of Grant et al.(2025)
                                                    # 1: run averaging on d_global_emergence to produce SI figure of Grant et al.(2025) of emergence fractions. Only activate when extr = 'all' (to verify later)

    flags['gdp_deprivation'] = 0                    # 0: do not process/load lifetime GDP/GRDI average
                                                    # 1: load lifetime GDP average analysis     
                                        
    flags['vulnerability'] = 0                      # 0: do not process subsets of d_collect_emergence vs gdp & deprivation quantiles
                                                    # 1: process/load d_collect_emergence vs gdp & deprivation quantiles for vulnerability analysis

    #----------------------------------------------------------#
    # Flags - Outputs                                          #
    #----------------------------------------------------------#

    flags['plots'] = 1                              # 0 do not produce and save plots 
                                                    # 1 produce and load plots 

    flags['reporting'] = 0                          # 0 do not produce results for reporting 
                                                    # 1 produce results for reporting
    
    # Use for specific climate extreme jobs - HPC only
    env_value_extr = os.getenv("EXTR_VALUE")
    if env_value_extr:  
        flags["extr"] = env_value_extr

        print(f"Using extr value: {flags['extr']}")

#%%----------------------------------------------------------------------------------------------------------#
# Execution - Initialisation and Execution of the subscripts based on the configuration                      #
#------------------------------------------------------------------------------------------------------------#

print("------------------------------------------------------------------------------------------------------")
print("|                            Start to run the Source2Suffering Project                               |")
print("------------------------------------------------------------------------------------------------------")
print("Current date and time:", datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
start_time = time.time() # Record the start time of execution
if env_value_paper:
    print("-------------------------------------------------------------------------------")
    print(f"|Model configuration - {env_value_paper}")
    print("-------------------------------------------------------------------------------")
else:
    print("-------------------------------------------------------------------------------")
    print(f"Model configuration - Manual")
    print("-------------------------------------------------------------------------------")
for key, value in flags.items():
    print(f"{key}: {value}\n")

print("-------------------------------------------------------------------------------")
print("Start to import settings")

from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

set_extremes(flags) # set extremes based on flag (this needs to happen here as it uses the flags dict defined above)

print("Settings imported")

print("-------------------------------------------------------------------------------")
print("Start to import and manipulate Demographic, GMT and ISIMIP data")

adr_load_manip = scripts_dir+"/load_manip.py"
with open(adr_load_manip) as f:
    exe_load_manip = f.read()
exec(exe_load_manip)

#%%------------------------------------------------------------------------------------
#                           Lifetime Exposure framework                               #
#--------------------------------------------------------------------------------------

print("-------------------------------------------------------------------------------")
print("Start Lifetime Exposure framework")

adr_exposure = scripts_dir+"/exposure.py"
with open(adr_exposure) as f:
    exe_exposure = f.read()
exec(exe_exposure)

#%%------------------------------------------------------------------------------------
#                       Emergence Lifetime Exposure framework                         #
#--------------------------------------------------------------------------------------

if Grant_2025==True:

    print("-------------------------------------------------------------------------------")
    print("Start Emergence Lifetime Exposure framework")

    adr_emergence = scripts_dir+"/emergence.py"
    with open(adr_emergence) as f:
        exe_emergence = f.read()
    exec(exe_emergence)

#%%------------------------------------------------------------------------------------
#                                Grid Scale Emergence                                 #
#--------------------------------------------------------------------------------------

if Grant_2025==True:

    print("-------------------------------------------------------------------------------")
    print("Start Gridscale Emergence framework")

    adr_gridscale = scripts_dir+"/gridscale.py"
    with open(adr_gridscale) as f:
        exe_gridscale = f.read()
    exec(exe_gridscale)


#%%------------------------------------------------------------------------------------
#                                   Outputs - Plots                                   #
#--------------------------------------------------------------------------------------

if flags['plots']:

    print("-------------------------------------------------------------------------------")
    print("Start Plots framework")
    print("-------------------------------------------------------------------------------")

    adr_plots = scripts_dir+"/plots.py"
    with open(adr_plots) as f:
        exe_plots = f.read()
    exec(exe_plots)

else : 
    print("No plots performed and saved")

#%%------------------------------------------------------------------------------------
#                                 Outputs - Reporting                                 #
#--------------------------------------------------------------------------------------

if flags['reporting']:


    print("-------------------------------------------------------------------------------")
    print("Start Reporting framework")
    print("-------------------------------------------------------------------------------")

    adr_report = scripts_dir+"/reporting.py"
    with open(adr_report) as f:
        exe_report = f.read()
    exec(exe_report)

else : 
    print("No reporting performed and saved")

#%%------------------------------------------------------------------------------------
#                                      Conclusion                                     #
#--------------------------------------------------------------------------------------
print("------------------------------------------------------------------------------------------------------")
print("|                        End of computations for the Source2Suffering Project                        |")
print("------------------------------------------------------------------------------------------------------")
# Calculate the script's execution time
end_time = time.time()
execution_time = end_time - start_time
# Convert execution time to hours
execution_time_in_hours = execution_time / 3600

# Print the execution time
print("Script execution time: {:.1f} seconds".format(execution_time))
print("Script execution time: {:.1f} minutes".format(execution_time/60))
print("Script execution time: {:.1f} hours".format(execution_time_in_hours))
print("------------------------------------------------------------------------------------------------------")
