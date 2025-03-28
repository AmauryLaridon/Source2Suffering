# --------------------------------------------------------------------------- #
# Subscript to execute the functions to load and manipulate data              #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------- #
# Libraries                                                       #
# --------------------------------------------------------------- #
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import pickle as pk
from scipy import interpolate
import regionmask
import glob
import os
from copy import deepcopy as cp

from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

sys.path.append(os.path.abspath(scripts_dir+"/pf_scripts"))
from pf_load_manip import *

# --------------------------------------------------------------- #
# Load global mean temperature projections                        #
# --------------------------------------------------------------- #


global df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_strj

df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_strj = load_GMT(
    year_start,
    year_end,
    year_range,
    flags,
)

print("GMT projections loaded")

# --------------------------------------------------------------- #
# Load and manipulate life expectancy, cohort and mortality data  #
# --------------------------------------------------------------- #

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

# --------------------------------------------------------------- #
# load ISIMIP model data                                          #
# --------------------------------------------------------------- #

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
ncountries = df_countries.shape[0] # number of available contries for the assessment
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


