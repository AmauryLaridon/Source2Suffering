# --------------------------------------------------------------------------- #
# Subscript to execute the functions to load and manipulate data              #
# --------------------------------------------------------------------------- #

#%%-------------------------------------------------------------- #
# Libraries                                                       #
# --------------------------------------------------------------- #

import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import regionmask as rm
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

if flags['gridscale_country_subset']:

    #--------------WORK IN PROGRESS for country_subset---------#

    # df_countries = d_countries['info_pop']
    # print("Under development gridscale_country_subset - Apply the framework for heatwavedarea analysis only on specific countries : {}".format(countries))                                                                   # In case we would like to perform the analysis only on a subset of coutries define
    # ind_countries = df_countries.reset_index(drop=True).index[df_countries["name"] == countries][0]      # in settings.py and pf_gridscale.py we need to retrieve the indices of that country for outputs
    # df_countries = df_countries.loc[df_countries['name'] == countries]

    # gdf_country_borders = d_countries['borders']
    # gdf_country_borders = gdf_country_borders.iloc[ind_countries]

    # da_population = d_countries['population_map']

    # df_birthyears = d_countries['birth_years']
    # print(type(df_birthyears))
    # print(df_birthyears.loc[ind_countries])
    # df_birthyears = df_birthyears.loc[ind_countries]

    # df_life_expectancy_5 = d_countries['life_expectancy_5']
    # df_life_expectancy_5 = df_life_expectancy_5[df_life_expectancy_5['name'] == countries]

    # da_cohort_size = d_countries['cohort_size']
    # da_cohort_size = da_cohort_size[ind_countries,:,:]

    # countries_regions, countries_mask = d_countries['mask']  

    #------------------------------------------------------------#

    df_countries = d_countries['info_pop']
    gdf_country_borders = d_countries['borders']
    da_regions = df_countries['region'].unique()
    da_population = d_countries['population_map']
    df_birthyears = d_countries['birth_years']
    df_life_expectancy_5 = d_countries['life_expectancy_5']
    da_cohort_size = d_countries['cohort_size']
    countries_regions, countries_mask = d_countries['mask']

else: 

    df_countries = d_countries['info_pop']
    gdf_country_borders = d_countries['borders']
    da_regions = df_countries['region'].unique()
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

global nruns, ncountries, nregions, nyears

nruns = len(d_isimip_meta)                      # number of available impact models runs used for this extreme
ncountries = df_countries.shape[0]              # number of available contries for the assessment  
nregions = len(da_regions)                      # number of regions 
nyears = len(year_range)                        # number of years for the assessment


# --------------------------------------------------------------- #
# load AR6 Regions                                                #
# --------------------------------------------------------------- #

grid_area = xr.open_dataarray(data_dir+'isimip/grid_resolution/clm45_area.nc4')

# arrays of lat/lon values
lat = grid_area.lat.values
lon = grid_area.lon.values

# 3d mask for ar6 regions. Gives a True of the gridcell fall into the specific region
da_ar6_regs_3D = rm.defined_regions.ar6.land.mask_3D(lon,lat)
# names of the ar6 regions 
da_ar6_regions_names = da_ar6_regs_3D.names.values
da_ar6_regions_index = np.arange(0,len(da_ar6_regions_names),1)
# create dictionnary with the index as key and names of the regions as values
d_ar6_regions = dict(zip(da_ar6_regions_index, da_ar6_regions_names))


# 3d mask for ar6 countries. Gives a True of the gridcell fall into the specific country
da_ar6_countries_3D = rm.mask_3D_geopandas(gdf_country_borders.reset_index(),lon,lat)
# indices of the 177 ar6 countries
da_ar6_countries_names = da_ar6_countries_3D.region.values

# --------------------------------------------------------------- #
# GMT steps tool                                                  #
# --------------------------------------------------------------- #

# stores for each GMT steps how many and which isimip simulations are available for remaping #
# only used for analysis 

sims_per_step = {}
for step in GMT_labels:
    sims_per_step[step] = []
    for i in list(d_isimip_meta.keys()):
        if d_isimip_meta[i]['GMT_strj_valid'][step]:
            sims_per_step[step].append(i)

print("ISMIP data loaded")


