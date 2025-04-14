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
import sys
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

    print('Processing country and regions data')

    d_countries ,df_regions, worldbank, unwpp = all_country_data(flags)

    print('Country and regions data loaded')

else: # load processed country data

    print('Loading processed country and regions data')

    # load country pickle
    d_countries = pk.load(open(data_dir+'{}/country/country_info.pkl'.format(flags['version']), 'rb'))

    # load regions pickle
    df_regions = pk.load(open(data_dir+'{}/country/regions_info.pkl'.format(flags['version']), 'rb'))

    # load worldbank pickle
    worldbank = pk.load(open(data_dir+'{}/country/worldbank.pkl'.format(flags['version']), 'rb'))

    # load unwpp pickle
    unwpp = pk.load(open(data_dir+'{}/country/unwpp.pkl'.format(flags['version']), 'rb'))

    print('Country and regions data loaded')
    
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

# print(type(da_cohort_size))
# print(np.shape(da_cohort_size))
# print(da_cohort_size)
# print("---------------")

d_cohort_size = get_cohortsize_countries(df_countries)

# print(type(d_cohort_size))
# print(np.shape(d_cohort_size))
# print(d_cohort_size["Zimbabwe"])

#sys.exit(0)
# --------------------------------------------------------------- #
# load Regions                                                    #
# --------------------------------------------------------------- #

#------------------------------- Manual configuration of ds_regions ------------------------------------#
# based on ms_manip.m from Thiery et al.(2021)                                                          #
#-------------------------------------------------------------------------------------------------------#

from scipy.io import loadmat

d_regions = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/regions_original.mat',squeeze_me=True)

d_regions = {k: v for k, v in d_regions.items() if not k.startswith('__')}

nregions = len(d_regions['name'])
coord_region = np.arange(0,nregions,1)

birth_year_val = np.arange(year_start,year_ref+1,1)
coord_birth_year = np.arange(0,year_ref+1-year_start,1)

ncountries = df_countries.shape[0]              # number of available contries for the assessment  
coord_country = np.arange(0,ncountries,1)

#------------------------------- Initialization of ds_regions -----------------------------------------#

# construction for the simple data variables that can be loaded from regions_original.mat
ds_regions = xr.Dataset(

            data_vars={
                'abbreviation' : (['region'], d_regions['abbreviation']),
                'name': (['region'], d_regions['name']),
                'birth_years': (['region','birth_year'], np.tile(birth_year_val, (nregions, 1))),
                'ind_member_countries': (['region','country'], np.zeros((nregions, ncountries), dtype=int)),
                #'member_countries', will be added to ds_regions() based on the definition of member_countries
            },

            coords={
                'region': coord_region,
                'birth_year': coord_birth_year,
                'country': coord_country,
            }

        )

# print("--------------")
# print(ds_regions)
# print("--------------")
# print(ds_regions["ind_member_countries"].values)
# print(ds_regions["ind_member_countries"].shape)
# print(ds_regions["ind_member_countries"].dtype)
# print("--------------")
# print(ds_regions["birth_years"].values)
# print(ds_regions["birth_years"].shape)
# print(ds_regions["birth_years"].dtype)

#d_regions_dvp = {} 


#-------- Configuration of ds_regions['ind_member_countries'] and ds_regions['member_countries'] ------#

member_countries = []
ind_member_countries = []

for i in range(len(coord_region)):

    # recover the name of the region in the loop over regions
    region_name = ds_regions['name'].loc[{'region': coord_region[i]}].values  # Assurez-vous de prendre la valeur correcte

    # condition of the region name is 'World' then recover all countries
    if region_name=='World':

        member_countries_all = df_countries['name'].tolist()
        ind_all = np.ones(ncountries, dtype=bool)

        member_countries.append(member_countries_all)
        ind_member_countries.append(ind_all)


    else:

        # find the countries which have an associated region equal to the one in the loop
        region_condition = df_countries['region'] == region_name
        income_condition = df_countries['incomegroup'] == region_name  
        
        # combine conditions for the name of the regions that could be based on geography or income
        matching_countries = df_countries[region_condition | income_condition].index.tolist()
    
        # add the name of the countries to the list
        member_countries.append(matching_countries)

        matching_index = df_countries[region_condition | income_condition].index
        
        # validation test that the list contains int
        if not np.issubdtype(matching_index.dtype, np.integer):
            matching_index = df_countries.index.get_indexer(matching_index)
        
        matching_index = matching_index.tolist()

        # create a boolean array of size ncountries initialized to False
        ind_array = np.zeros(ncountries, dtype=bool)
        ind_array[matching_index] = True

        # add to the final list
        ind_member_countries.append(ind_array)

# conversation to array of the lists
member_countries = np.array(member_countries, dtype=object)
ind_member_countries = np.array(ind_member_countries, dtype=object)

# print(member_countries)
# print(type(member_countries))
# print(np.shape(member_countries))

# print("------------")

# print(ind_member_countries)
# print(type(ind_member_countries))
# print(np.shape(ind_member_countries))

#--------- Integration into ds_regions['ind_member_countries'] and ds_regions['member_countries']------#

for i in range(nregions):
    # Affectation du masque booléen des pays membres pour la région i
    ds_regions['ind_member_countries'][i, :] = ind_member_countries[i].astype(int)

ds_regions['member_countries'] = xr.DataArray(
    data=member_countries,
    dims=['region'],
    coords={'region': coord_region}
)


# print("--------------")
# print(ds_regions)
# print("--------------")
# print(ds_regions["member_countries"][0].values)
# print(ds_regions["member_countries"][0].shape)
# print(ds_regions["member_countries"][0].dtype)

# print("--------------")
# print(ds_regions["member_countries"][1].values)
# print(ds_regions["member_countries"][1].shape)
# print(ds_regions["member_countries"][1].dtype)
# print("--------------")
# print(ds_regions["member_countries"].values)
# print(ds_regions["member_countries"].shape)
# print(ds_regions["member_countries"].dtype)
# print(ds_regions["member_countries"][0].values)

#--------------------------- Configuration of ds_regions['cohort_weights']  ---------------------------#

# sys.exit(0)

# tmp1 = []
# tmp2 = []
# tmp1 = 0

# da_cohort_size[country,ind_2020,ages]


# Liste pour accumuler les cohort_weights par région
# all_weights = []
# region_labels = []

# for region_name in regions["name"]:
#     # Récupérer les pays membres de la région
#     ind_member_countries = regions["members"][region_name]

#     tmp_list = []
#     for country in ind_member_countries:
#         da = countries["cohort_size"][country]  # DataArray: [time, age]
#         # Extraire la cohorte née en year_ref
#         da_ref = da.sel(time=year_ref)
#         tmp_list.append(da_ref)

#     # Combiner en un seul DataArray: (country, age)
#     tmp2 = xr.concat(tmp_list, dim="country")
#     tmp2 = tmp2.assign_coords(country=ind_member_countries)

#     # Réordonner pour avoir (age, country)
#     cohort_w = tmp2.transpose("age", "country")

#     # Ajouter une nouvelle dimension région
#     cohort_w = cohort_w.expand_dims(region=[region_name])
#     all_weights.append(cohort_w)
#     region_labels.append(region_name)

# Concaténer toutes les régions sur une seule dimension "region"
# cohort_weights = xr.concat(all_weights, dim="region")
# cohort_weights = cohort_weights.assign_coords(region=region_labels)


#------------------------------- Luke's get_regions_data() -----------------------------#


df_worldbank_region = worldbank[1]
df_unwpp_region = unwpp[1]

d_region_countries, df_birthyears_regions, df_life_expectancy_5_regions, d_cohort_weights_regions = get_regions_data(
    df_countries, 
    df_regions, 
    df_worldbank_region, 
    df_unwpp_region, 
    d_cohort_size,
)

# print(type(d_cohort_weights_regions['South Asia']))
# print(np.shape(d_cohort_weights_regions['South Asia']))
# print("d_cohort_weights_regions", d_cohort_weights_regions['South Asia'])

#sys.exit(0)

# ---------------------------------- Intégration dans ds_regions['cohort_weights'] ---------------------------------- #

# --- Option will all countries as coord --- #

# # Liste pour stocker les DataArrays de chaque région
# da_list = []

# # Liste complète des pays comme dans ds_regions
# all_countries = ds_regions.coords['country'].values

# # On boucle sur chaque région selon l'ordre des index dans ds_regions
# for i in range(len(ds_regions.coords['region'])):

#     region_name = ds_regions['name'].sel(region=i).item()
#     print(f"Ajout des poids pour la région : {region_name}")

#     df_region = d_cohort_weights_regions[region_name].copy()
#     df_region.columns.name = None  # enlever le nom des colonnes

#     # Réindexer les colonnes pour inclure tous les pays
#     df_region_full = df_region.reindex(columns=all_countries)

#     # Définir les coordonnées 'age'
#     age_coords = df_region_full.index

#     # Création du DataArray
#     da = xr.DataArray(
#         data=df_region_full.values,
#         dims=["age", "country"],
#         coords={
#             "age": age_coords,
#             "country": df_region_full.columns,
#             "region": i
#         },
#         name="cohort_weights"
#     )

#     da_list.append(da)

# # Concaténation sur la dimension région
# cohort_weights_all = xr.concat(da_list, dim="region")

# # Ajout au dataset principal
# ds_regions["cohort_weights"] = cohort_weights_all

# print("--------------")
# print(ds_regions)
# print("--------------")
# print(ds_regions['cohort_weights'].sel(region=7,country=25))
# # cohort_weights_na = ds_regions['cohort_weights'].sel(region=7)
# # print(cohort_weights_na)
# # print(cohort_weights_na.shape)
# # print(type(cohort_weights_na))
# print("--------------")


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

global nruns, nyears

nruns = len(d_isimip_meta)                      # number of available impact models runs used for this extreme
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


