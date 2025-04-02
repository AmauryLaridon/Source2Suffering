# --------------------------------------------------------------------------- #
# Subscript to execute the gridscale analysis                                 #
# --------------------------------------------------------------------------- #
            
#%%-------------------------------------------------------------- #
# Libraries                                                       #
# --------------------------------------------------------------- #

import os
import glob
import requests
from zipfile import ZipFile
import io
import xarray as xr
import pickle as pk
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
import matplotlib as mpl
import mapclassify as mc
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import regionmask as rm
import numpy as np
import pandas as pd
import geopandas as gpd
# import rioxarray as rxr
from scipy import interpolate
from scipy.stats import ttest_rel
import cartopy.crs as ccrs
from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

sys.path.append(os.path.abspath(scripts_dir+"/pf_scripts"))
from pf_gridscale import *

# --------------------------------------------------------------- #
# Execute gridscale                                               #
# --------------------------------------------------------------- #

# list of countries to run gridscale analysis on (sometimes doing subsets across basiss/regions in floods/droughts)
gridscale_countries = get_gridscale_regions(
    grid_area,
    flags,
    gdf_country_borders,
)

# print(gridscale_countries)
# print(type(gridscale_countries))
# print(np.shape(gridscale_countries))

# birth year aligned cohort sizes for gridscale analysis (summed over lat/lon per country)
if not os.path.isfile(data_dir+'{}/country/gs_cohort_sizes.pkl'.format(flags['version'])):

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
