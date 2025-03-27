#!/usr/bin/env python3
# -----------------------------------------------------------------------------------------------------------
# Main script to postprocess and visualise lifetime exposure to climate extremes data
#
# Python translation of the MATLAB scripts of Thiery et al. (2021)
# The translation of the scripts has been first performed by Luke Grant for Grant et al.(2025)
# The modification of Luke Grant's scripts have been performed by Amaury Laridon for Laridon et al.(2025)
# -----------------------------------------------------------------------------------------------------------
# Associated papers
# -----------------------------------------------------------------------------------------------------------
#
# - Thiery et al.(2021) - https://www.science.org/doi/abs/10.1126/science.abi7339
# - Grant et al.(2025) - in review
# - Laridon et al.(2025) - in prep
#
# -----------------------------------------------------------------------------------------------------------
# Summary and notes
# -----------------------------------------------------------------------------------------------------------
#
# Data types are defined in the variable names starting with:  
#     df_     : DataFrame    (pandas)
#     gdf_    : GeoDataFrame (geopandas)
#     da_     : DataArray    (xarray)
#     d_      : dictionary  
#     sf_     : shapefile
#     ...dir  : directory
# -----------------------------------------------------------------------------------------------------------
     
#%%------------------------------------------------------------------------------------
# Libraries
#--------------------------------------------------------------------------------------

import xarray as xr
import pickle as pk
import time
#from mpl_toolkits.axes_grid1 import make_axes_locatable
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

#%%------------------------------------------------------------------------------------
# Flags - Define the Configuration of the Framework 
#--------------------------------------------------------------------------------------

# Global flags use to configure the framework

global flags

flags = {}

# Boolean to produce the exact configuration of the framework and outputs to 
# reproduce the results of the associated paper

global Thiery_2021, Grant_2025, Laridon_2025

Thiery_2021 = False 
Grant_2025 = False      # define to True if you want to run parts of the framework that relies on the emergence and gridscale
Laridon_2025 = True
env_value_paper= 0 

#--------------------------------------------------------------------------------------
# Submit jobs that execute main.py with specific configurations - use on HPC only     #
#--------------------------------------------------------------------------------------

# Use for specific paper's configuration jobs
env_value_paper = os.getenv("CONFIG_PAPER_VALUE")
if env_value_paper:
    Thiery_2021 = False 
    Grant_2025 = False 
    Laridon_2025 = False 
    if env_value_paper == "thiery_2021":
        Thiery_2021 = True
    if env_value_paper == "grant_2025":
        Grant_2025 = True
    if env_value_paper == "laridon_2025":
        Laridon_2025 = True    

#%%------------------------------------------------------------------------------------
# Configuration of the Framework to reproduce the papers
#--------------------------------------------------------------------------------------

if Thiery_2021==True:
    flags['extr'] = 'heatwavedarea'
    flags['gmt'] = 'original'
    flags['rm'] = 'no_rm'
    flags['version'] = 'pickles'
    flags['run'] = 1
    flags['mask'] = 1
    flags['lifetime_exposure_cohort'] = 1
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
    flags['extr'] = 'heatwavedarea'
    flags['gmt'] = 'ar6_new'
    flags['rm'] = 'rm'
    flags['version'] = 'pickles_v3'
    flags['run'] = 0
    flags['mask'] = 0
    flags['lifetime_exposure_cohort'] = 0
    flags['lifetime_exposure_pic'] = 0
    flags['emergence'] = 0
    flags['birthyear_emergence'] = 0
    flags['gridscale'] = 1
    flags['gridscale_le_test'] = 0
    flags['gridscale_country_subset'] = 0
    flags['global_emergence_recollect'] = 0
    flags['global_avg_emergence'] = 0
    flags['gdp_deprivation'] = 0
    flags['vulnerability'] = 0
    flags['plots'] = 1
    flags['reporting'] = 1

if Laridon_2025==True:
    print('Configuration of the Framework for Laridon et al.(2025) not settle going to Manual Configuration')

#%%------------------------------------------------------------------------------------
# Flags - Manual Configuration of the Framework
#--------------------------------------------------------------------------------------

if not env_value_paper:

    flags['extr'] = 'heatwavedarea'                # 0: all
                                                    # 1: burntarea
                                                    # 2: cropfailedarea
                                                    # 3: driedarea
                                                    # 4: floodedarea
                                                    # 5: heatwavedarea
                                                    # 6: tropicalcyclonedarea

    flags['gmt'] = 'ar6_new'                            # original: use Wim's stylized trajectory approach with max trajectory a linear increase to 3.5 deg                               
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

    flags['lifetime_exposure_cohort'] = 0           # 0: do not process ISIMIP runs to compute exposure across cohorts (i.e. load exposure pickle)
                                                    # 1: process ISIMIP runs to compute exposure across cohorts (i.e. produce and save exposure as pickle)   
                                                                        
    flags['lifetime_exposure_pic'] = 1              # 0: do not process ISIMIP runs to compute picontrol exposure (i.e. load exposure pickle)
                                                    # 1: process ISIMIP runs to compute picontrol exposure (i.e. produce and save exposure as pickle)

    flags['emergence'] = 0                          # 0: do not process ISIMIP runs to compute cohort emergence (i.e. load cohort exposure pickle)
                                                    # 1: process ISIMIP runs to compute cohort emergence (i.e. produce and save exposure as pickle)

    #--------------------- Produce Error ----------------------

    flags['birthyear_emergence'] = 0                # 0: only run calc_birthyear_align with birth years from 1960-2020
                                                    # 1: run calc_birthyear_align with birth years from 1960-2100. Produces an error.  
                                                
    #----------------------------------------------------------

                            
    flags['gridscale'] = 0                          # 0: do not process grid scale analysis, load pickles
                                                    # 1: process grid scale analysis

    flags['gridscale_le_test'] = 0                  # 0: do not process the grid scale analysis testing diff versions of constant life expectancy
                                                    # 1: process grid scale analysis testing diff versions of constant life expectancy    
                                                                
    flags['gridscale_country_subset'] = 0           # 0: run gridscale analysis on all countries
                                                    # 1: run gridscale analysis on subset of countries determined in "get_gridscale_regions" 

    #--------------------- Produce Error ----------------------

    flags['global_emergence_recollect'] = 0         # 0: do not load pickles of global emergence masks used for vulnerability assessment
                                                    # 1: load pickles                  

    flags['global_avg_emergence'] = 0               # 0: do not run d_global_emergence used for SI Figures of Grant et al.(2025)
                                                    # 1: run averaging on d_global_emergence to produce SI figure of Grant et al.(2025) of emergence fractions. Only activate when extr = 'all' (to verify later)

    flags['gdp_deprivation'] = 0                    # 0: do not process/load lifetime GDP/GRDI average
                                                    # 1: load lifetime GDP average analysis     
                                        
    flags['vulnerability'] = 0                      # 0: do not process subsets of d_collect_emergence vs gdp & deprivation quantiles
                                                    # 1: process/load d_collect_emergence vs gdp & deprivation quantiles for vulnerability analysis

    #----------------------------------------------------------
    # Flags - Outputs
    #----------------------------------------------------------

    flags['plots'] = 0                               # 0 do not produce and save plots 
                                                     # 1 produce and load plots 

    flags['reporting'] = 1                          # 0 do not produce results for reporting 
                                                    # 1 produce results for reporting
    
    # Use for specific climate extreme jobs - HPC only
    env_value_extr = os.getenv("EXTR_VALUE")
    if env_value_extr:  
        flags["extr"] = env_value_extr

        print(f"Using extr value: {flags['extr']}")

#--------------------------------------------------------------------------------------
# Init
#--------------------------------------------------------------------------------------

print("-----------------------------------------------------------")
print("         Start to run the Source2Suffering Project")
print("-----------------------------------------------------------")
print("Current date and time:", datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
start_time = time.time() # Record the start time of execution
if env_value_paper:
    print("-----------------------------------------------------------")
    print(f"Model configuration - {env_value_paper}")
    print("-----------------------------------------------------------")
else:
    print("-----------------------------------------------------------")
    print(f"Model configuration - Manual")
    print("-----------------------------------------------------------")
for key, value in flags.items():
    print(f"{key}: {value}\n")
#%%------------------------------------------------------------------------------------
# Settings
#--------------------------------------------------------------------------------------
print("--------------------------------------------------")
print("Start to import settings")
print("--------------------------------------------------")


from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins = init()

# set extremes based on flag (this needs to happen here as it uses the flags dict defined above)
set_extremes(flags)

print("Settings imported")

#%%------------------------------------------------------------------------------------
# Load and manipulate demographic, GMT and ISIMIP data
#--------------------------------------------------------------------------------------
print("--------------------------------------------------")
print("Start to import Demographic, GMT and ISIMIP data")
print("--------------------------------------------------")


from load_manip import *

# ------------------------------------------------------------------
# Load global mean temperature projections
# ------------------------------------------------------------------

global df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_strj

df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_strj = load_GMT(
    year_start,
    year_end,
    year_range,
    flags,
)

print("GMT projections loaded")

# ------------------------------------------------------------------
# Load and manipulate life expectancy, cohort and mortality data
# ------------------------------------------------------------------

if flags['mask']: # compute and process country info

    print('Processing country data')

    d_countries = all_country_data(flags)

    print('Country data loaded')

else: # load processed country data

    print('Loading processed country data')

    # load country pickle
    d_countries = pk.load(open(data_dir+'{}/country/country_info.pkl'.format(flags['version']), 'rb'))

    print('Country data loaded')
    
# unpack country information
df_countries = d_countries['info_pop']
gdf_country_borders = d_countries['borders']
da_population = d_countries['population_map']
df_birthyears = d_countries['birth_years']
df_life_expectancy_5 = d_countries['life_expectancy_5']
da_cohort_size = d_countries['cohort_size']
countries_regions, countries_mask = d_countries['mask']  

# ------------------------------------------------------------------
# load ISIMIP model data
# ------------------------------------------------------------------

d_isimip_meta,d_pic_meta = load_isimip(
    extremes,
    model_names,
    df_GMT_15,
    df_GMT_20,
    df_GMT_NDC,
    df_GMT_strj,
    flags,
)

global nruns, ncountries, nyears

nruns = len(d_isimip_meta) # number of available impact models runs used for this extreme
ncountries = np.shape(df_countries[0]) # number of available contries for the assessment
nyears = len(year_range) # number of years for the assessment

# stores for each GMT steps how many and which isimip simulations are available for remaping #
# only used for analysis 

sims_per_step = {}
for step in GMT_labels:
    sims_per_step[step] = []
    for i in list(d_isimip_meta.keys()):
        if d_isimip_meta[i]['GMT_strj_valid'][step]:
            sims_per_step[step].append(i)

print("ISMIP data loaded")


#%%------------------------------------------------------------------------------------
# Lifetime Exposure framework
#--------------------------------------------------------------------------------------

print("--------------------------------------------------")
print("Start Lifetime Exposure framework")
print("--------------------------------------------------")

from exposure import *

# ------------------------------------------------------------------
# process lifetime exposure across cohorts
# ------------------------------------------------------------------

if flags['lifetime_exposure_cohort']:

    print("Computing Lifetime Exposure across cohorts")
    start_time = time.time()
    
    
    # calculate exposure per country and per cohort
    calc_cohort_lifetime_exposure(
        d_isimip_meta,
        df_countries,
        countries_regions,
        countries_mask,
        da_population,
        da_cohort_size,
        flags,
    )
    
    print("--- {} minutes to compute Lifetime Exposure across Cohorts for all countries ---".format(
        np.floor((time.time() - start_time) / 60),
        )
          )

# ------------------------------------------------------------------
# load lifetime exposure across cohorts
# ------------------------------------------------------------------

    
else: # load processed cohort exposure data
    
    print('Loading processed Lifetime Exposures across cohorts. !Not yet settle!')

    # for i in len(d_isimip_meta):

          # add a way to store those but needs to be uptaded later

    #     with open(data_dir+'{}/{}/exposure_peryear_perage_percountry_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
    #         d_exposure_perrun_pic = pk.load(f)

# --------------------------------------------------------------------
# process lifetime exposure across cohorts for PIC climate conditions
# --------------------------------------------------------------------

if flags['lifetime_exposure_pic']:
    
    print('Computing PIC Lifetime Exposures across cohorts')
    start_time = time.time()
    
    d_exposure_perrun_pic = calc_lifetime_exposure_pic(
        d_pic_meta, 
        df_countries, 
        countries_regions, 
        countries_mask, 
        da_population, 
        df_life_expectancy_5, 
        flags,
    )
    
    print("--- {} minutes to compute Lifetime Exposure across Cohorts for all countries under PIC climate conditions for the 1960 demographics ---".format(
        np.floor((time.time() - start_time) / 60),
        )
          )    

# --------------------------------------------------------------------
# load lifetime exposure across cohorts for PIC climate conditions
# --------------------------------------------------------------------

else:  
    
    print('Loading PIC Lifetime Exposures across cohorts')

    with open(data_dir+'{}/{}/exposure_pic_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
        d_exposure_perrun_pic = pk.load(f)

# --------------------------------------------------------------------
# Multi-model mean of exposure
# --------------------------------------------------------------------

# computation of multi-model mean of the exposure
ds_exposure_pic = calc_exposure_mmm_pic_xr(
    d_exposure_perrun_pic,
    'country',
    'pic',
)

#%%------------------------------------------------------------------------------------
# Emergence Lifetime Exposure framework
#--------------------------------------------------------------------------------------

if Grant_2025==True:

    print("--------------------------------------------------")
    print("Start Emergence Lifetime Exposure framework")
    print("--------------------------------------------------")


    from emergence import *

    global grid_area
    grid_area = xr.open_dataarray(data_dir+'isimip/clm45_area.nc4')

    if flags['emergence']:

        # --------------------------------------------------------------------
        # process emergence of cumulative exposures, mask cohort exposures for time steps of emergence

        print("Computing Emergence of cumulative exposures")
        
        if flags['birthyear_emergence']:
            
            by_emergence = np.arange(1960,2101)
            
        else:
            
            by_emergence = birth_years        
        
        if not os.path.isfile(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version'])):
            
            # need new cohort dataset that has total population per birth year (using life expectancy info; each country has a different end point)
            da_cohort_aligned = calc_birthyear_align(
                da_cohort_size,
                df_life_expectancy_5,
                by_emergence,
            )
            
            # convert to dataset and add weights
            ds_cohorts = ds_cohort_align(
                da_cohort_size,
                da_cohort_aligned,
            )
            
            # pickle birth year aligned cohort sizes and global mean life expectancy
            with open(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version']), 'wb') as f:
                pk.dump(ds_cohorts,f)  

        else:
            
            # load pickled birth year aligned cohort sizes and global mean life expectancy
            with open(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version']), 'rb') as f:
                ds_cohorts = pk.load(f)                             
        
        ds_ae_strj, ds_pf_strj = strj_emergence(
            d_isimip_meta,
            df_life_expectancy_5,
            ds_exposure_pic,
            ds_cohorts,
            by_emergence,
            flags,
        )
            
    else: # load pickles

        print("Loading Emergence of cumulative exposures")
        
        pass
        
        # birth year aligned population
        with open(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version']), 'rb') as f:
            ds_cohorts = pk.load(f)
        
        # pop frac
        with open(data_dir+'{}/{}/pop_frac_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
            ds_pf_strj = pk.load(f)              
    

#%%------------------------------------------------------------------------------------
# Grid Scale Emergence
#--------------------------------------------------------------------------------------

if Grant_2025==True:

    print("--------------------------------------------------")
    print("Start Gridscale Emergence framework")
    print("--------------------------------------------------")

    from gridscale import *

    # list of countries to run gridscale analysis on (sometimes doing subsets across basiss/regions in floods/droughts)
    gridscale_countries = get_gridscale_regions(
        grid_area,
        flags,
        gdf_country_borders,
    )

    # birth year aligned cohort sizes for gridscale analysis (summed over lat/lon per country)
    if not os.path.isfile(data_dir+'{}/country/gs_cohort_sizes.pkl'.format(flags['version'])):

        #print('getting da_gs_popdenom')
        da_gs_popdenom = get_gridscale_popdenom(
            gridscale_countries,
            da_cohort_size,
            countries_mask,
            countries_regions,
            da_population,
            df_life_expectancy_5,
        )

        # pickle birth year aligned cohort sizes for gridscale analysis (summed per country)
        with open(data_dir+'{}/country/gs_cohort_sizes.pkl'.format(flags['version']), 'wb') as f:
            pk.dump(da_gs_popdenom,f)  
            
    else:
        
        # load pickle birth year aligned cohort sizes for gridscale analysis (summed per country, i.e. not lat/lon explicit)
        #print('loading da_gs_popdenom')
        with open(data_dir+'{}/country/gs_cohort_sizes.pkl'.format(flags['version']), 'rb') as f:
            da_gs_popdenom = pk.load(f)

    # run gridscale emergence analysis
    if flags['gridscale']:
                
        print("Computing Gridscale Emergence of cumulative exposures")
        
        ds_pf_gs = gridscale_emergence(
            d_isimip_meta,
            d_pic_meta,
            flags,
            gridscale_countries,
            da_cohort_size,
            countries_regions,
            countries_mask,
            df_life_expectancy_5,
            da_population,
        )  
        
    else:
        
        print("Loading Gridscale Emergence of cumulative exposures")
        # # load pickled aggregated lifetime exposure, age emergence and pop frac datasets
        # with open(data_dir+'{}/{}/gridscale_aggregated_lifetime_exposure_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
        #     ds_le_gs = pk.load(f)
        with open(data_dir+'{}/{}/gridscale_aggregated_pop_frac_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
            ds_pf_gs = pk.load(f)

            
    if flags['gridscale_le_test']:

        print("Computing Test of Gridscale Emergence of cumulative exposures")
        
        ds_pf_gs_le_test = gridscale_emergence_life_expectancy_constant(
            d_isimip_meta,
            d_pic_meta,
            flags,
            gridscale_countries,
            da_cohort_size,
            countries_regions,
            countries_mask,
            df_life_expectancy_5,
            da_population,
        )        
        
    else:
        
        print("Loading Test of Gridscale Emergence of cumulative exposures")

        with open(data_dir+'{}/{}/gridscale_aggregated_pop_frac_le_test_{}.pkl'.format(flags['version'],flags['extr']+'_le_test',flags['extr']), 'rb') as f:
            ds_pf_gs_le_test = pk.load(f)    
                

    # read in all global emergence masks (d_global_emergence is then used for vulnerability assessment, but only possible on hpc because it is large for some hazards)
    if flags['global_emergence_recollect']:

        print("--------------------------------------------------")
        print("Start Global Emergence Mask for Vulnerability framework")
        print("--------------------------------------------------")

        from gridscale import *

        # temporarily commented out extremes in this function outside heatwaved area to test new means extraction below
        d_global_emergence = collect_global_emergence(
            grid_area,
            flags,
            countries_mask,
            countries_regions,
            gridscale_countries,
            df_GMT_strj,
        )
        
        # temporarily commented out extremes in this function outside heatwaved area to test new means extraction below
        d_global_pic_qntls = collect_pic_qntls(
            grid_area,
            flags,
            gridscale_countries,
            countries_mask,
            countries_regions,
        )  
        
        d_global_pic_qntls_extra = collect_pic_qntls_extra(
            grid_area,
            flags,
            gridscale_countries,
            countries_mask,
            countries_regions,
        )    
        

    if flags['global_avg_emergence']:

        print("--------------------------------------------------")
        print("Start Averaging of Emergence framework")
        print("--------------------------------------------------")

        from gridscale import *
        
        # run averaging on d_global_emergence to produce SI figure of emergence fractions
        ds_emergence_mean = get_mean_emergence(
            df_GMT_strj,
            flags,
            da_population,
            d_global_emergence,
        )    
        
    # load/proc GDP and deprivation data
    if flags['gdp_deprivation']:

        print("--------------------------------------------------")
        print("Start GDP/GRDI framework")
        print("--------------------------------------------------")

        from gridscale import *
        
        ds_gdp, ds_grdi = load_gdp_deprivation(
            flags,
            grid_area,
            da_population,
            countries_mask,
            countries_regions,
            gridscale_countries,
            df_life_expectancy_5,
        )
        
    # vulnerability subsetting
    if flags['vulnerability']:  

        print("--------------------------------------------------")
        print("Start Vulnerability framework")
        print("--------------------------------------------------")

        from gridscale import *

        # get spatially explicit cohort sizes for all birth years in analysis
        da_cohort_size_1960_2020 = get_spatially_explicit_cohorts_1960_2020(
            flags,
            gridscale_countries,
            countries_mask,
            countries_regions,
            da_cohort_size,
            da_population,
        )
        
        # adds data arrays to ds_gdp and ds_grdi with ranked vulnerability binned by population (i.e. ranges of ranked vulnerability, physically distributed, grouped/binned by population size)            
        ds_gdp_qntls, ds_grdi_qntls = get_vulnerability_quantiles(
            flags,
            grid_area,
            da_cohort_size_1960_2020,
            ds_gdp,
            ds_grdi,
        )
            
        # just a dummy d_global_emergence to run emergence_by_vulnerability
        # try:
        #     d_global_emergence
        # except NameError:
        #     print('to save memory on my laptop, d_global_emergence is not unpickled. defining a dummy var for emergence_by_vulnerability')
        #     d_global_emergence={}
        # else:
        #     pass

        # dataset of emergence numbers selected by quantiles of vulnerability, both with grdi and gdp
        ds_vulnerability = emergence_by_vulnerability(
            flags,
            df_GMT_strj,
            ds_gdp_qntls,
            ds_grdi_qntls,
            da_cohort_size_1960_2020,
            d_global_emergence,
        )
    

#%%------------------------------------------------------------------------------------
# Outputs - Plots
#--------------------------------------------------------------------------------------

if flags['plots']:

    print("--------------------------------------------------")
    print("Start Plots framework")
    print("--------------------------------------------------")

    adr_plots = scripts_dir+"/plots.py"
    with open(adr_plots) as f:
        exe_plots = f.read()
    exec(exe_plots)

else : 
    print("No plots performed and saved")

#%%------------------------------------------------------------------------------------
# Outputs - Reporting
#--------------------------------------------------------------------------------------

if flags['reporting']:


    print("--------------------------------------------------")
    print("Start Reporting framework")
    print("--------------------------------------------------")

    adr_report = scripts_dir+"/reporting.py"
    with open(adr_report) as f:
        exe_report = f.read()
    exec(exe_report)

else : 
    print("No reporting performed and saved")

#%%------------------------------------------------------------------------------------
# Conclusion
#--------------------------------------------------------------------------------------

print("-----------------------------------------------------------")
print("   End of computations for the Source2Suffering Project")
print("-----------------------------------------------------------")
# Calculate the script's execution time
end_time = time.time()
execution_time = end_time - start_time
# Convert execution time to hours
execution_time_in_hours = execution_time / 3600

# Print the execution time
print("Script execution time: {:.1f} seconds".format(execution_time))
print("Script execution time: {:.1f} minutes".format(execution_time/60))
print("Script execution time: {:.1f} hours".format(execution_time_in_hours))
print("-----------------------------------------------------------")
