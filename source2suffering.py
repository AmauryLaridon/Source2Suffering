# ----------------------------------------------------------------------------------------
# Script that defines the framework for the Source2Suffering project. This defines the 
# functions neede to assess the additionnal number of people facing one additational
# climate extreme across their lifetime due to a certain amount of CO2 emissions.
# Also compute heat-related mortality associated with the emissions using the 
# concept of 'mortality cost of carbon' from Bessler et al.(2021) in Nature Com.
# ----------------------------------------------------------------------------------------
#%%  ----------------------------------------------------------------
# Libraries  
# -------------------------------------------------------------------

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
from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()


# -------------------------------------------------------------------
# Init
# -------------------------------------------------------------------


# Define Transient Climate response to cumulative emissions # 
# Expert advice (27/03/2024): You can use 1.65°C per 1000 PgC, but (if useful) you can also use transparently communicated higher percentiles.
# These provide a risk perspective. For example, 2.2°C per 1000 PgC is still only the 66th percentile, which is not necessarily extreme.
# However, this method is mainly applicable to CO₂ and not CO₂-eq. So if there is a lot of methane in those emissions, you would need to adjust it slightly.

TCRE = 0.45/1000E9 # 1.65C / 3.7 = 0,45C per 1000 Gt CO2eq


# Define the mortality cost of Carbon
mortality_cost_carbon = 4434 #4,434 metric tons of carbon dioxide in 2020 [] causes one excess death globally in expectation between 2020-2100 (Bressler, 2021)


#%%------------------------------------------------------------------------------------
# Transcription of the mf_emissions2npeople.m script from Wim Thiery to python. This
# function will serve as a basis for the Source2Suffering development below. 
# 
# Framework to compute additional number of people faceing one additational
# climate extreme across their lifetime due to a certain amount of CO2 emissions.
# Also compute heat-related mortality associated with the emissions using the 
# concept of 'mortality cost of carbon' from Bessler et al.(2021) in Nature Com. 
#--------------------------------------------------------------------------------------


def emissions2npeople(CO2_emissions, TCRE, exposure_perregion_BE, birth_years, year_start, year_end, GMT_BE, valp_cohort_size_abs, rounding, ind_extreme):
    """Compute the number of people affected by additional climate extremes in their lifetime
    and the number of heat-related deaths between today and 2100."""

    # Compute change in GMT from emissions 
    dGMT = TCRE * CO2_emissions                           

    # Generate list of birth years for iteration
    years_loop = list(range(year_end, year_start - 1, -1))
    nbirthyears = len(years_loop)

    # Initialize arrays
    nr_newborns = np.zeros(nbirthyears)
    nr_extra_climate_extremes_newborns = np.zeros(nbirthyears)
    nr_children_facing_extra_climate_extreme = np.zeros(nbirthyears)

    # Loop over birth years from year_end to year_start
    for i in range(nbirthyears):

        # Extract lifetime exposure data and GMT anomaly in 2100 #
        idx = np.where(birth_years == years_loop[i])[0]
        valc_exposure_climate_extreme_newborns = np.squeeze(exposure_perregion_BE[ind_extreme, 12, idx, :]) # lifetime absolute climate extreme exposure for particular generation
        valc_GMT_2100 = np.squeeze(GMT_BE[-1, :])  # 2100 GMT anomaly

        # Fit a linear curve between the exposure to climate extreme and the GMT anomaly and extract the slope #
        valc_pf = np.polyfit(valc_GMT_2100, valc_exposure_climate_extreme_newborns, 1)
        valc_slope_exposure_climate_extreme = valc_pf[0]

        # Extract the number of people in the cohort #
        nr_newborns[i] = valp_cohort_size_abs[-i-1, 12]

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

    # Compute total for 2015-2020 #
    nr_children_facing_extra_climate_extreme = np.append(
        nr_children_facing_extra_climate_extreme, 
        np.nansum(nr_children_facing_extra_climate_extreme)
    )

    return nr_children_facing_extra_climate_extreme