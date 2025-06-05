# --------------------------------------------------------------------------- #
# Subscript to execute the functions to load and manipulate data              #
# --------------------------------------------------------------------------- #

#%%---------------------------------------------------------------#
# Libraries                                                       #
# ----------------------------------------------------------------#

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

# --------------------------------------------------------------- #
# Execution of the sub_script                                     #
# --------------------------------------------------------------- #

adr_pf_exposure = scripts_dir+'/pf_scripts/pf_exposure.py'
with open(adr_pf_exposure) as f:
    exec(f.read(), globals())

# ----------------------------------------------------------------#
# Process Land Fraction Exposed (LFE) across ISIMIP runs          #
# ----------------------------------------------------------------#

if flags['landfraction_exposed']:

    print("----------------------------------------------------------------------")
    print("Computing Land Fraction Exposed (LFE) per ISIMIP run to {}       ".format(flags['extr']))
    print("---------------------------------------------------------------------")
    start_time = time.time()

    #--------------  Compute per ISIMIP run Land Fraction Exposed --------------#

    ds_lfe_percountry_perrun, ds_lfe_perregion_perrun = calc_landfraction_exposed(
        d_isimip_meta,
        df_countries,
        countries_regions,
        countries_mask,
        ds_regions,
        flags,)

    #-------------------  Compute MMM Land Fraction Exposed --------------------#

    print("\n---------------------------------------------------------------")
    print("Computing MMM Land Fraction Exposed (LFE) to {}              ".format(flags['extr']))
    print("---------------------------------------------------------------")

    ds_lfe_percountry = calc_landfraction_exposed_mmm(ds_lfe_percountry_perrun, flags)
    
    ds_lfe_perregion = calc_landfraction_exposed_mmm(ds_lfe_perregion_perrun, flags)

    print("\n--- {} minutes to compute Land Fraction Exposed for all countries and regions ---".format(
        np.floor((time.time() - start_time) / 60),
        )
          )

# ----------------------------------------------------------------#
# Load Land Fraction Exposed (LFE) across ISIMIP runs             #
# ----------------------------------------------------------------#
  
else: # load processed land fraction exposed data

    #---------------  Load per ISIMIP run Land Fraction Exposed -------------#

    print('Loading processed Land Fraction Exposed per ISIMIP simulation')

    with open(data_dir+'{}/{}/ds_lfe_percountry_perrun_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
        ds_lfe_percountry_perrun = pk.load(f)

    with open(data_dir+'{}/{}/ds_lfe_perregion_perrun_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
        ds_lfe_perregion_perrun = pk.load(f)

    #--------------------  Load MMM Land Fraction Exposed -------------------#

    print('\nLoading processed MMM Land Fraction Exposed\n')

    with open(data_dir+'{}/{}/ds_lfe_percountry_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
        ds_lfe_percountry = pk.load(f)

    with open(data_dir+'{}/{}/ds_lfe_perregion_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
        ds_lfe_perregion = pk.load(f)

# ----------------------------------------------------------------#
# Process lifetime exposure across cohorts                        #
# ----------------------------------------------------------------#

if flags['lifetime_exposure']:

    print("---------------------------------------------------------------")
    print("Computing Lifetime Exposure (LE) per ISIMIP run to {}       ".format(flags['extr']))
    print("---------------------------------------------------------------")
    start_time = time.time()

    if Thiery_2021 or Source2Suffering:

        #-----------------  Compute per ISIMIP run Lifetime Exposure ---------------#
        
        ds_le_percountry_perrun, ds_le_perregion_perrun = calc_lifetime_exposure(
            d_isimip_meta,
            df_countries,
            countries_regions,
            countries_mask,
            da_population,
            df_life_expectancy_5,
            ds_regions,
            da_cohort_size_regions,
            flags,)

        #---------------------  Compute MMM Lifetime Exposure ---------------------#

        print("\n---------------------------------------------------------------")
        print("Computing MMM Lifetime Exposure (LE) to {}              ".format(flags['extr']))
        print("---------------------------------------------------------------")

        ds_le_percountry = calc_lifetime_exposure_mmm(ds_le_percountry_perrun, flags)
        
        ds_le_perregion = calc_lifetime_exposure_mmm(ds_le_perregion_perrun, flags)

        #----------------------------- Compute EMF  -------------------------------#

        print("\n---------------------------------------------------------------")
        print("Computing EMF of Lifetime Exposure (LE) to {}            ".format(flags['extr']))
        print("---------------------------------------------------------------")

        ds_EMF_percountry = calc_EMF(flags, ds_le_exposure=ds_le_percountry, ref_pic = False)

        ds_EMF_perregion = calc_EMF(flags, ds_le_exposure=ds_le_perregion, ref_pic = False)

    
    if Grant_2025:
        
        print("Lifetime Exposure computations will be performed in the emergence analysis")

        # function developped by L.Grant for Grant et al.(2025) #
        # calculate exposure per country and per cohort to try to analyse the "age of emergence"
        # by being time/age explicit to assess this. Did not come to fruition and not retain for further
        # usage in Grant et al.(2025)

        calc_cohort_lifetime_exposure(
            d_isimip_meta,
            df_countries,
            countries_regions,
            countries_mask,
            da_population,
            da_cohort_size,
            flags,
        )
    
    print("\n--- {} minutes to compute Lifetime Exposure for all countries and regions ---".format(
        np.floor((time.time() - start_time) / 60),
        )
          )

# ----------------------------------------------------------------#
# Load lifetime exposure across cohorts                           #
# ----------------------------------------------------------------#

    
else: # load processed cohort exposure data

    if Grant_2025:
    
        print('Loading processed Lifetime Exposures across cohorts will be done in the emergence analysis')

    elif Thiery_2021 or Source2Suffering:
        
        #-----------------  Load per ISIMIP run Lifetime Exposure ---------------#

        print('Loading processed Lifetime Exposure per ISIMIP simulation')

        with open(data_dir+'{}/{}/ds_le_percountry_perrun_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_percountry_perrun = pk.load(f)
        
        with open(data_dir+'{}/{}/ds_le_perregion_perrun_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_perregion_perrun = pk.load(f)

        #---------------------  Load MMM Lifetime Exposure ---------------------#

        print('\nLoading processed MMM Lifetime Exposure')

        with open(data_dir+'{}/{}/ds_le_percountry_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_percountry = pk.load(f)
        
        with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_perregion = pk.load(f)

        #----------------------------- Load EMF  -------------------------------#

        print('\nLoading EMF of Lifetime Exposure')

        with open(data_dir+'{}/{}/ds_EMF_percountry_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:

            ds_EMF_percountry_GMT = pk.load(f)

        with open(data_dir+'{}/{}/ds_EMF_perregion_gmt_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt'],flags['rm']), 'rb') as f:

            ds_EMF_perregion_GMT = pk.load(f)

# --------------------------------------------------------------- #
# Process lifetime exposure across cohorts for                    #
# PIC climate conditions                                          #
# --------------------------------------------------------------- #

if flags['lifetime_exposure_pic']:

    #-----------------  Compute per ISIMIP run Lifetime Exposure ---------------#
    
    print("\n ------------------------------------------------------------------")
    print("|            Computing Lifetime Exposure to {} under PI         ".format(flags['extr']))
    print(" ------------------------------------------------------------------")
    start_time = time.time()
    
    d_pic_le_percountry_perrun = calc_lifetime_exposure_pic(
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
    
    print('\nLoading PIC Lifetime Exposure per ISIMIP simulation')

    with open(data_dir+'{}/{}/d_pic_le_percountry_perrun.pkl'.format(flags['version'],flags['extr']), 'rb') as f:
            d_pic_le_percountry_perrun = pk.load(f)

# --------------------------------------------------------------- #
# Multi-model mean of exposure for PIC simulations                #
# --------------------------------------------------------------- #

if Grant_2025:

    # computation of multi-model mean of the exposure
    ds_exposure_pic = calc_exposure_mmm_pic_xr(
        d_pic_le_percountry_perrun,
        'country',
        'pic',
    )
