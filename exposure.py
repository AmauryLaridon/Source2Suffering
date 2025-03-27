# --------------------------------------------------------------------------- #
# Subscript to execute the functions to load and manipulate data              #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------- #
# Libraries                                                       #
# --------------------------------------------------------------- #


import sys
from operator import index
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import pickle as pk
from scipy import interpolate
from scipy import stats as sts
import regionmask as rm
import glob
import time
import matplotlib.pyplot as plt
from copy import deepcopy as cp
from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins = init()

sys.path.append(os.path.abspath(scripts_dir+"/pf_scripts"))
from pf_exposure import *

# --------------------------------------------------------------- #
# Execute lifetime exposure across cohorts                        #
# --------------------------------------------------------------- #

if flags['lifetime_exposure_cohort']:

    print("Computing Lifetime Exposure across cohorts")
    start_time = time.time()
    

    if Grant_2025:
    
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

# --------------------------------------------------------------- #
# Load lifetime exposure across cohorts                           #
# --------------------------------------------------------------- #

    
else: # load processed cohort exposure data

    if Grant_2025:
    
        print('Loading processed Lifetime Exposures across cohorts. !Not yet settle!')

    else:

        print('Loading processed Lifetime Exposures across cohorts. !Not yet settle!')

        # for i in len(d_isimip_meta):

            # add a way to store those but needs to be uptaded later

        #     with open(data_dir+'{}/{}/exposure_peryear_perage_percountry_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
        #         d_exposure_perrun_pic = pk.load(f)

# --------------------------------------------------------------- #
# Process lifetime exposure across cohorts for                    #
# PIC climate conditions                                          #
# --------------------------------------------------------------- #

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

# --------------------------------------------------------------- #
# Load lifetime exposure across cohorts for                       #
# PIC climate conditions                                          #
# --------------------------------------------------------------- #

else:  
    
    print('Loading PIC Lifetime Exposures across cohorts')

    with open(data_dir+'{}/{}/exposure_pic_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
        d_exposure_perrun_pic = pk.load(f)

# --------------------------------------------------------------- #
# Multi-model mean of exposure                                    #
# --------------------------------------------------------------- #

if Grant_2025:
    # computation of multi-model mean of the exposure
    ds_exposure_pic = calc_exposure_mmm_pic_xr(
        d_exposure_perrun_pic,
        'country',
        'pic',
    )
