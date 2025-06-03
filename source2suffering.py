# ---------------------------------------------------------------------------------------- #
# Script that defines the framework for the Source2Suffering project. This defines the     #
# functions neede to assess the additionnal number of people facing one additational       #
# climate extreme across their lifetime due to a certain amount of CO2 emissions.          #
# Also compute heat-related mortality associated with the emissions using the              #
# concept of 'mortality cost of carbon' from Bessler et al.(2021) in Nature Com.           #
# ---------------------------------------------------------------------------------------- #
#%%  ---------------------------------------------------------------- #
# Libraries                                                           #
# ------------------------------------------------------------------- #

import os
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
import openpyxl
import pickle as pk
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import regionmask as rm
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import interpolate
import cartopy.crs as ccrs

#%%------------------------------------------------------------------------------------ #
# Transcription of the mf_emissions2npeople.m script from Wim Thiery to python. This    #
# function will serve as a basis for the Source2Suffering development below.            #
#                                                                                       #
# Framework to compute additional number of people faceing one additational             #
# climate extreme across their lifetime due to a certain amount of CO2 emissions.       #
# Also compute heat-related mortality associated with the emissions using the           #
# concept of 'mortality cost of carbon' from Bessler et al.(2021) in Nature Com.        #
#-------------------------------------------------------------------------------------- #


def emissions2npeople(CO2_emissions, TCRE, ds_le, region_ind, year_start, year_end, df_GMT_strj, da_valp_cohort_size_abs, rounding):
    """Compute the number of people affected by additional climate extremes in their lifetime
    due to specific CO2 emissions"""

    # Compute change in GMT from emissions 
    dGMT = TCRE * CO2_emissions                           

    # Generate list of birth years for iteration
    years_loop = list(range(year_end, year_start - 1, -1))
    nbirthyears = len(years_loop)

    # Initialize arrays
    nr_newborns = np.zeros(nbirthyears)
    nr_extra_climate_extremes_newborns = np.zeros(nbirthyears)
    nr_children_facing_extra_climate_extreme = np.zeros(nbirthyears)
    slope_exposure = np.zeros(nbirthyears)

    # Loop over birth years from year_end to year_start
    for i in range(nbirthyears):

        # Extract lifetime exposure data and GMT anomaly in 2100 
        # for each cohort among the birth_years between year_start and year_end #
        valc_exposure_climate_extreme_newborns = ds_le['mmm_BE'].sel(region=region_ind, birth_year=years_loop[i])

        birth_year=years_loop[i]

        # Recover the 2113 GMT anomaly ? To be confirmed. Should be placed outside of the loop
        valc_GMT_2100 = df_GMT_strj.iloc[-1]  

        # print(valc_exposure_climate_extreme_newborns)
        # print(np.shape(valc_exposure_climate_extreme_newborns))

        # print(valc_GMT_2100)
        # print(np.shape(valc_GMT_2100))

        # Fit a linear curve between the exposure to climate extreme and the GMT anomaly and extract the slope #
        valc_pf = np.polyfit(valc_GMT_2100, valc_exposure_climate_extreme_newborns, 1)
        valc_slope_exposure_climate_extreme = valc_pf[0]
        slope_exposure[i] = valc_slope_exposure_climate_extreme

        # Extract the number of people in the cohort #

        if year_end == 2020:

            #nr_newborns[i] = da_valp_cohort_size_abs[-(1+i), region_ind]
            nr_newborns[i] = da_valp_cohort_size_abs.sel(region=region_ind, ages=ages[-(1+i)]).item()

            #print("Number of people in 2020 of the birth years = {} - nr_newborns[i] = {}".format(years_loop[i],nr_newborns[i]))
        
        if year_start == 1960:

            #nr_newborns[i] = da_valp_cohort_size_abs[(year_end-year_start)-i, region_ind]
            nr_newborns[i] = da_valp_cohort_size_abs.sel(region=region_ind, ages=ages[(year_end-year_start)-i]).item()

            #print("Number of people in 2020 of the birth years = {} - nr_newborns[i] = {}".format(years_loop[i],nr_newborns[i]))

        # Compute the average change in lifetime exposure #
        nr_extra_climate_extremes_newborns[i] = valc_slope_exposure_climate_extreme * dGMT

        # Compute the number of children facing at least one additional extreme #
        val = nr_newborns[i] * nr_extra_climate_extremes_newborns[i]

        if rounding == 0:
            nr_children_facing_extra_climate_extreme[i] = np.floor(val)
        elif rounding == 1:
            nr_children_facing_extra_climate_extreme[i] = np.floor(val / 1000) * 1000
        elif rounding == 2:
            nr_children_facing_extra_climate_extreme[i] = np.floor(val / 100) * 100

        # Ensure no negative values #
        nr_children_facing_extra_climate_extreme[i] = max(nr_children_facing_extra_climate_extreme[i], 0)

    # Compute total for all the birth cohorts of interest #
    nr_children_facing_extra_climate_extreme = np.append(
        nr_children_facing_extra_climate_extreme, 
        np.nansum(nr_children_facing_extra_climate_extreme)
    )

    return nr_children_facing_extra_climate_extreme, slope_exposure


#%%------------------------------------------------------------------------------------ #
# Compute heat-related mortality associated with the emissions using the                #
# concept of 'mortality cost of carbon' from Bessler et al.(2021) in Nature Com.        #
#-------------------------------------------------------------------------------------- #

def mortality_cost_carbon(CO2_emissions):

    # 4,434 metric tons of carbon dioxide in 2020 cfr. (Bressler, 2021) 
    # causes one excess death globally in expectation between 2020-2100 

    mortality_cost_carbon_val = 4434 

    valc_mortality = (CO2_emissions / mortality_cost_carbon_val // 1000) * 1000

    return valc_mortality