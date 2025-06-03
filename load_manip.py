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

# --------------------------------------------------------------- #
# Execution of the sub_script                                     #
# --------------------------------------------------------------- #

adr_pf_load_manip = scripts_dir+'/pf_scripts/pf_load_manip.py'
with open(adr_pf_load_manip) as f:
    exec(f.read(), globals())

# --------------------------------------------------------------- #
# Load global mean temperature projections                        #
# --------------------------------------------------------------- #


global df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_OS, df_GMT_noOS, ds_GMT_STS, df_GMT_strj 

df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_OS, df_GMT_noOS, ds_GMT_STS, df_GMT_strj = load_GMT(
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

d_cohort_size = get_cohortsize_countries(df_countries,flags)

# Save as pickle
with open(data_dir + '{}/country/df_countries.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(df_countries, f)
with open(data_dir + '{}/country/gdf_country_borders.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(gdf_country_borders, f)
with open(data_dir + '{}/country/da_regions.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(da_regions, f)
with open(data_dir + '{}/country/da_population.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(da_population, f)
with open(data_dir + '{}/country/df_birthyears.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(df_birthyears, f)
with open(data_dir + '{}/country/df_life_expectancy_5.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(df_life_expectancy_5, f)
with open(data_dir + '{}/country/da_cohort_size.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(da_cohort_size, f)
with open(data_dir + '{}/country/countries_regions.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(countries_regions, f)
with open(data_dir + '{}/country/countries_mask.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(countries_mask, f)
with open(data_dir + '{}/country/d_cohort_size.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(d_cohort_size, f)

# --------------------------------------------------------------------------------------- #
# Build ds_regions partialy based on loading the original object from Thiery et al.(2021) #
# based on ms.manip.m                                                                     #
# --------------------------------------------------------------------------------------- #

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
                'mask': (['region'], d_regions['mask'])
            },

            coords={
                'region': coord_region,
                'birth_year': coord_birth_year,
                'country': coord_country,
            }

        )

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

#--------- Integration into ds_regions['ind_member_countries'] and ds_regions['member_countries']------#

for i in range(nregions):
    # Affectation du masque booléen des pays membres pour la région i
    ds_regions['ind_member_countries'][i, :] = ind_member_countries[i].astype(int)

ds_regions['member_countries'] = xr.DataArray(
    data=member_countries,
    dims=['region'],
    coords={'region': coord_region}
)

with open(data_dir + '{}/country/ds_regions.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(ds_regions, f)

# --------------------------------------------------------------------- #
# Construction of d_cohort_weights_regions based on                     #
# d_cohort_size and the get_regions_data() function written by          #
# L.Grant but that was not used in final analysis of Grant et al.(2025) #
# --------------------------------------------------------------------- #

df_worldbank_region = worldbank[1]
df_unwpp_region = unwpp[1]

d_region_countries, df_birthyears_regions, df_life_expectancy_5_regions, d_cohort_weights_regions = get_regions_data(
    df_countries, 
    df_regions, 
    df_worldbank_region, 
    df_unwpp_region, 
    d_cohort_size,
    flags,
)

# --------------------------------------------------------------------- #
# Construction of da_cohort_size_regions based on da_cohort_size which  #
# is used in the final analysis of Grant et al.(2025).                  #
# The get_regions_cohort() function is written by A.Laridon for         # 
# Laridon et al.(2025)                                                  #
# --------------------------------------------------------------------- #

da_cohort_size_regions = get_regions_cohort(
    df_countries, 
    ds_regions, 
    da_cohort_size, 
    flags
)


# --------------------------------------------------------------------- #
# Construction of valp_cohort_size_abs which containts the total        # 
# number of population of a given age per region. This object is        # 
# used in source2suffering.py                                           #
# --------------------------------------------------------------------- #

if Source2Suffering:

    # --------------------------------------------------------------- #
    # Construction of valp_cohort_size_abs based on ms_valp.m from    #
    # Thiery et al.(2021).                                            #
    # --------------------------------------------------------------- #
    
    regions_list = list(d_cohort_weights_regions.keys())
    nregions = len(regions_list)

    valp_cohort_size_abs = np.zeros((len(ages), nregions))

    # Loop on the ages 
    for ind_age, age in enumerate(ages):
        
        # Loop over regions
        for ind_region, region_name in enumerate(regions_list):
            
            # Access the DataFrame for the region
            df_cohort = d_cohort_weights_regions[region_name]  # shape: [age x countries]
            
            # Check that the age is in the index
            if age in df_cohort.index:
                # Select the row corresponding to the age (series with countries as columns)
                cohort_row = df_cohort.loc[age]
                
                # Sum across countries and convert to absolute values
                total = cohort_row.sum() * 1000
                
                # Store the rounded value
                valp_cohort_size_abs[ind_age, ind_region] = np.round(total, 1)
            else:
                # If age is missing, set as NaN
                valp_cohort_size_abs[ind_age, ind_region] = np.nan
    
if Thiery_2021:

    # --------------------------------------------------------------- #
    # Importation of valp_cohort_size_abs from Thiery et al.(2021)    #                                            #
    # --------------------------------------------------------------- #

    valp_cohort_size_abs = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/valp_cohort_size_abs.mat',squeeze_me=True)
    valp_cohort_size_abs = valp_cohort_size_abs['valp_cohort_size_abs']

# -------------------------------------- #
# Save the object as DataSet in a pickle #
# -------------------------------------- #

regions = np.arange(nregions)  # regions indexed from 0 to nregions-1

# Create a DataArray
da_valp_cohort_size_abs = xr.DataArray(
    data=valp_cohort_size_abs,  # shape (age, region)
    dims=("age", "region"),
    coords={
        "age": ages,
        "region": regions
    },
    name="cohort_size_abs"
)

# Create the Dataset
ds_valp_cohort_size_abs = xr.Dataset({"cohort_size_abs": da_valp_cohort_size_abs})

# Save as pickle
with open(data_dir + '{}/country/ds_valp_cohort_size_abs.pkl'.format(flags['version']), 'wb') as f:
    pk.dump(ds_valp_cohort_size_abs, f)

# --------------------------------------------------------------- #
# load ISIMIP model data                                          #
# --------------------------------------------------------------- #

d_isimip_meta,d_pic_meta = load_isimip(
    extremes,
    model_names,
    df_GMT_15,
    df_GMT_20,
    df_GMT_NDC,
    df_GMT_OS,
    df_GMT_noOS,
    ds_GMT_STS,
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


