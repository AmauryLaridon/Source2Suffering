# ------------------------------------------------------------------------------ #
# Functions needed and developped for the Source2Suffering project               #
#                                                                                #
# This script contains in its first part the technical functions such            #
# as emissions2npeople(), monte_carlo_sampling() and the functions to builds     #
# PDFs for the 4 main variables needed to link emissions to human                #
# lifetime exposure.                                                             #
#                                                                                #
# In te second part, this script contains the functions developped for specific  #
# computations for Laridon et al.(2025) and for all the Assessment Reports       # 
# based on the Source2Suffering Framework.                                       #
# ------------------------------------------------------------------------------ #

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


#-------------------------------------------------------------------------------------------------------------------------------#
#                                           Source2Suffering Technical Core Functions                                           #
#-------------------------------------------------------------------------------------------------------------------------------#


#%%------------------------------------------------------------------------------------ #
# Transcription of the mf_emissions2npeople.m script from Wim Thiery to python. This    #
# function will serve as a basis for the Source2Suffering development below.            #
#                                                                                       #
# Framework to compute additional number of people faceing one additational             #
# climate extreme across their lifetime due to a certain amount of CO2 emissions.       #
#-------------------------------------------------------------------------------------- #


def emissions2npeople(CO2_emissions, TCRE, ds_le, region_ind, year_start, year_end, df_GMT_strj, da_valp_cohort_size_abs, rounding):
    """Compute the number of people affected by additional climate extremes in their lifetime
    due to specific CO2 emissions"""

    # Compute change in GMT from emissions 
    dGMT = TCRE * CO2_emissions             

    # List of the ages of interest
    ages = np.arange(60,-1,-1)              

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

            nr_newborns[i] = da_valp_cohort_size_abs.sel(region=region_ind, ages=ages[-(1+i)]).item()
        
        if year_start == 1960:

            nr_newborns[i] = da_valp_cohort_size_abs.sel(region=region_ind, ages=ages[(year_end-year_start)-i]).item()

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
        elif rounding == 3:
            nr_children_facing_extra_climate_extreme[i] = val

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


#-------------------------------------------------------------------------------------------------------------------------------#
#                                   Source2Suffering Applications for Laridon et al.(2025)                                      #
#-------------------------------------------------------------------------------------------------------------------------------#

def reference_pulse():

    print("\n ---------------------------------------------------------")
    print("|          Assessment for a Reference Pulse               |")
    print("|                      of 1 GtCO2                         |")
    print(" ---------------------------------------------------------")

    # -------------------------------------------------------------------------- #
    # Define Total emissions of the reference fossil fuel project under study    #
    # -------------------------------------------------------------------------- #
    
    # TOTAL emissions reference pulse (MtCO2e: add E6 to express as tCO2e) #
    CO2_emissions_Reference_Pulse = 1000e6

    # -------------------------------------------------------------------------- #
    # Define Transient Climate response to cumulative emission                   #
    # -------------------------------------------------------------------------- #

    # Use of the reference value given by the inital expert advice (27/03/24)

    TCRE_init = 0.45 / 1e12 # 1.65°C / 3.7 = 0,45°C per 1000 Gt CO2eq
    
    # -------------------------------------------------------------------------- #
    # GMT produce by Reference Pulse using the TCRE                              #
    # -------------------------------------------------------------------------- #

    print("\nCO2 emissions of the Reference Pulse = {} MtCO2eq\n".format(CO2_emissions_Reference_Pulse/10**6))

    print("Value used for TCRE = 0.45 °C per 1000 Gt CO2eq\n")

    print("GMT produce by Reference Pulse = {} °C".format(TCRE_init*CO2_emissions_Reference_Pulse))

    # -------------------------------------------------------------------------- #
    # Baseline Question : How many additional hazards will the 2010-2020         #
    # birth cohort experience due to the emissions of the Reference Pulse        # 
    # compared to the reference birth cohorts of 1960-1970                       #
    # -------------------------------------------------------------------------- #

    print("")
    print(" --------------------------------------------------------------------------------------------------------------------")
    print("| Q.(1) - Number of children in the world facing an additional hazard due to the total emissions the Reference Pulse |")
    print(" --------------------------------------------------------------------------------------------------------------------")
    print("")

    year_start_as = 2010
    year_end_as = 2020
    year_start_as_ref = 1960
    year_end_as_ref = 1970

    all_extremes = [
            'burntarea', 
            'cropfailedarea', 
            'driedarea', 
            'floodedarea', 
            'heatwavedarea', 
            'tropicalcyclonedarea'
        ]

    hazards_name = [
        'Wild Fires', 
        'Crop Failures', 
        'Droughts', 
        'River Floods', 
        'Heatwaves', 
        'Tropical Cyclones'
    ]

    # Dictionary to store values for each extreme hazard type
    d_valc = {}
    d_valc_ref = {}
    d_valc_total = {}
    d_valc_total_ref = {}
    d_valc_slope_exposure = {}
    d_valc_slope_exposure_ref = {}

    n = 0
    # Loop over all extreme hazard types
    for extr in all_extremes:

        extr_name = hazards_name[n]
        print("---------- Hazard = {} ----------\n".format(extr_name))
        n+=1

        # Load the corresponding exposure dataset with new demography
        with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_perregion = pk.load(f)

        # Load the absolute cohort sizes at the regional level with the new demography
        with open(data_dir + '{}/country/da_valp_cohort_size_abs.pkl'.format(flags['version']), 'rb') as f:
            da_valp_cohort_size_abs = pk.load(f)

        # Compute the number of children exposed to the extra hazard under the Reference Pulse
        valc_nr_children_facing_extra_hazard_Reference_Pulse, S2S_slope_exposure = emissions2npeople(
            CO2_emissions=CO2_emissions_Reference_Pulse,
            TCRE=TCRE_init,
            ds_le=ds_le_perregion,
            region_ind=11,
            year_start=year_start_as,
            year_end=year_end_as,
            df_GMT_strj=df_GMT_strj,
            da_valp_cohort_size_abs=da_valp_cohort_size_abs,
            rounding=2
        )

        # Compute the number of children exposed to the extra hazard under the Reference Pulse scenario for the reference birth cohorts
        valc_nr_children_facing_extra_hazard_Reference_Pulse_ref, S2S_slope_exposure_ref = emissions2npeople(
            CO2_emissions=CO2_emissions_Reference_Pulse,
            TCRE=TCRE_init,
            ds_le=ds_le_perregion,
            region_ind=11,
            year_start=year_start_as_ref,
            year_end=year_end_as_ref,
            df_GMT_strj=df_GMT_strj,
            da_valp_cohort_size_abs=da_valp_cohort_size_abs,
            rounding=2
        )

        nr_total = valc_nr_children_facing_extra_hazard_Reference_Pulse[-1]
        nr_total_ref = valc_nr_children_facing_extra_hazard_Reference_Pulse_ref[-1]

        # Store the value of the total number of people exposed
        d_valc_total[extr] = nr_total
        d_valc_total_ref[extr] = nr_total_ref

        # Store the value of slope exposure for each hazard
        d_valc_slope_exposure[extr] = S2S_slope_exposure
        d_valc_slope_exposure_ref[extr] = S2S_slope_exposure_ref

        # Exclude the last value along the only dimension (since the array is 1D)
        da_trimmed = valc_nr_children_facing_extra_hazard_Reference_Pulse[:-1]
        da_ref_trimmed = valc_nr_children_facing_extra_hazard_Reference_Pulse_ref[:-1]

        # Reverse the values
        da_reversed = da_trimmed[::-1]
        da_ref_reversed = da_ref_trimmed[::-1]

        # Store the result in the dictionary under the key 'extr'
        d_valc[extr] = da_reversed
        d_valc_ref[extr] = da_ref_reversed

        # Generate list of birth years in descending order
        years_loop = list(range(year_end_as, year_start_as - 1, -1))
        nbirthyears = len(years_loop)
        years_loop_ref = list(range(year_end_as_ref, year_start_as_ref - 1, -1))
        nbirthyears_ref = len(years_loop_ref)

        # Print number of exposed children per birth year
        for i in range(nbirthyears):
            print("For birth year {} = {} children".format(
                years_loop[i], int(valc_nr_children_facing_extra_hazard_Reference_Pulse[i])))

        # Print total number of exposed children across all birth years
        print("For total birth of the {}–{} period = {} children \n".format(
            year_start_as, year_end_as, int(valc_nr_children_facing_extra_hazard_Reference_Pulse[-1]))
        )
        # Print total number of exposed children across all birth years for the reference period
        print("For total birth of the {}–{} period = {} people \n".format(
            year_start_as_ref, year_end_as_ref, int(valc_nr_children_facing_extra_hazard_Reference_Pulse_ref[-1]))
        )

    # ----------------------- Creating DataArray and DataSet ---------------------------- #

    # Create a DataArray containing all values for each hazard and birth year
    da_valc_nr_children_facing_extra_hazard_Reference_Pulse = xr.DataArray(
        data=[d_valc[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as, year_end_as + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_nr_children_facing_extra_hazard_Reference_Pulse"
    )

    da_valc_nr_children_facing_extra_hazard_Reference_Pulse_ref = xr.DataArray(
        data=[d_valc_ref[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as_ref, year_end_as_ref + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_nr_children_facing_extra_hazard_Reference_Pulse_ref"
    )

    # Create a DataArray containing the values of the total of people exposed for each hazard
    da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse = xr.DataArray(
        data=[d_valc_total[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
        },
        dims=["hazard"],
        name="valc_total_nr_children_facing_extra_hazard_Reference_Pulse"
    )

    da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref = xr.DataArray(
        data=[d_valc_total_ref[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
        },
        dims=["hazard"],
        name="valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref"
    )

    # Create a DataArray containing the values of the slope exposure for each hazard

    da_valc_slope_exposure = xr.DataArray(
        data=[d_valc_slope_exposure[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as, year_end_as + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_slope_exposure"
    )

    da_valc_slope_exposure_ref = xr.DataArray(
        data=[d_valc_slope_exposure_ref[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as_ref, year_end_as_ref + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_slope_exposure_ref"
    )

    ds_S2S_Reference_Pulse = xr.Dataset(
    {
        "valc_nr_children_facing_extra_hazard_Reference_Pulse": da_valc_nr_children_facing_extra_hazard_Reference_Pulse,
        "valc_total_nr_children_facing_extra_hazard_Reference_Pulse": da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse,
        "valc_slope_exposure": da_valc_slope_exposure,
        "valc_nr_children_facing_extra_hazard_Reference_Pulse_ref": da_valc_nr_children_facing_extra_hazard_Reference_Pulse_ref,
        "valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref": da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref,
        "valc_slope_exposure_ref": da_valc_slope_exposure_ref

    }
    )

    # -------------------------------------------------------------------------- #
    # Spatialization Question : How many does the exposure spread among space ?  #
    # -------------------------------------------------------------------------- #

    print("")
    print(" ---------------------------------------------------- ")
    print("| Q.(2) - How does the exposure spread among space ? |")
    print(" ---------------------------------------------------- ")
    print("")

    # Define the list of geographical regions
    region_int_ind = [0, 1, 3, 6, 7, 8, 9, 11] 
    region_names = ds_regions['name'].sel(region=region_int_ind).values

    # Prepare empty lists to store datasets for each region
    list_da_valc_nr = []
    list_da_valc_nr_ref = []
    list_da_valc_total = []
    list_da_valc_total_ref = []
    list_da_valc_slope = []
    list_da_valc_slope_ref = []

    for region_idx, region_name in zip(region_int_ind, region_names):
        print(f"\n------------------------------------------------")
        print(f"   Computing for REGION = {region_name}")
        print(f"------------------------------------------------\n")

        d_valc = {}
        d_valc_ref = {}
        d_valc_total = {}
        d_valc_total_ref = {}
        d_valc_slope_exposure = {}
        d_valc_slope_exposure_ref = {}

        n = 0
        for extr in all_extremes:
            extr_name = hazards_name[n]
            print(f"Hazard = {extr_name}\n")
            n += 1

            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_perregion = pk.load(f)

            with open(data_dir + '{}/country/da_valp_cohort_size_abs.pkl'.format(flags['version']), 'rb') as f:
                da_valp_cohort_size_abs = pk.load(f)

            # Exposures for current birth cohort
            valc_nr, slope_expo = emissions2npeople(
                CO2_emissions=CO2_emissions_Reference_Pulse,
                TCRE=TCRE_init,
                ds_le=ds_le_perregion,
                region_ind=region_idx,
                year_start=year_start_as,
                year_end=year_end_as,
                df_GMT_strj=df_GMT_strj,
                da_valp_cohort_size_abs=da_valp_cohort_size_abs,
                rounding=3
            )

            # Exposures for reference birth cohort
            valc_nr_ref, slope_expo_ref = emissions2npeople(
                CO2_emissions=CO2_emissions_Reference_Pulse,
                TCRE=TCRE_init,
                ds_le=ds_le_perregion,
                region_ind=region_idx,
                year_start=year_start_as_ref,
                year_end=year_end_as_ref,
                df_GMT_strj=df_GMT_strj,
                da_valp_cohort_size_abs=da_valp_cohort_size_abs,
                rounding=3
            )

            d_valc_total[extr] = valc_nr[-1]
            d_valc_total_ref[extr] = valc_nr_ref[-1]
            d_valc_slope_exposure[extr] = slope_expo
            d_valc_slope_exposure_ref[extr] = slope_expo_ref
            d_valc[extr] = valc_nr[:-1][::-1]
            d_valc_ref[extr] = valc_nr_ref[:-1][::-1]

        # Create DataArrays for this region
        da_valc_nr = xr.DataArray(
            data=[d_valc[extr] for extr in all_extremes],
            coords={"hazard": all_extremes, "birth_year": list(range(year_start_as, year_end_as + 1))},
            dims=["hazard", "birth_year"]
        )

        da_valc_nr_ref = xr.DataArray(
            data=[d_valc_ref[extr] for extr in all_extremes],
            coords={"hazard": all_extremes, "birth_year": list(range(year_start_as_ref, year_end_as_ref + 1))},
            dims=["hazard", "birth_year"]
        )

        da_valc_total = xr.DataArray(
            data=[d_valc_total[extr] for extr in all_extremes],
            coords={"hazard": all_extremes},
            dims=["hazard"]
        )

        da_valc_total_ref = xr.DataArray(
            data=[d_valc_total_ref[extr] for extr in all_extremes],
            coords={"hazard": all_extremes},
            dims=["hazard"]
        )

        da_slope = xr.DataArray(
            data=[d_valc_slope_exposure[extr] for extr in all_extremes],
            coords={"hazard": all_extremes, "birth_year": list(range(year_start_as, year_end_as + 1))},
            dims=["hazard", "birth_year"]
        )

        da_slope_ref = xr.DataArray(
            data=[d_valc_slope_exposure_ref[extr] for extr in all_extremes],
            coords={"hazard": all_extremes, "birth_year": list(range(year_start_as_ref, year_end_as_ref + 1))},
            dims=["hazard", "birth_year"]
        )

        list_da_valc_nr.append(da_valc_nr)
        list_da_valc_nr_ref.append(da_valc_nr_ref)
        list_da_valc_total.append(da_valc_total)
        list_da_valc_total_ref.append(da_valc_total_ref)
        list_da_valc_slope.append(da_slope)
        list_da_valc_slope_ref.append(da_slope_ref)

    # Convert lists into stacked DataArrays using region names as coordinate
    da_valc_nr_children_facing_extra_hazard_Reference_Pulse = xr.concat(list_da_valc_nr, dim="region")
    da_valc_nr_children_facing_extra_hazard_Reference_Pulse["region"] = region_names

    da_valc_nr_children_facing_extra_hazard_Reference_Pulse_ref = xr.concat(list_da_valc_nr_ref, dim="region")
    da_valc_nr_children_facing_extra_hazard_Reference_Pulse_ref["region"] = region_names

    da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse = xr.concat(list_da_valc_total, dim="region")
    da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse["region"] = region_names

    da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref = xr.concat(list_da_valc_total_ref, dim="region")
    da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref["region"] = region_names

    da_valc_slope_exposure = xr.concat(list_da_valc_slope, dim="region")
    da_valc_slope_exposure["region"] = region_names

    da_valc_slope_exposure_ref = xr.concat(list_da_valc_slope_ref, dim="region")
    da_valc_slope_exposure_ref["region"] = region_names

    # Assemble the final dataset with named regions
    ds_S2S_Reference_Pulse_Regions = xr.Dataset(
        {
            "valc_nr_children_facing_extra_hazard_Reference_Pulse": da_valc_nr_children_facing_extra_hazard_Reference_Pulse,
            "valc_nr_children_facing_extra_hazard_Reference_Pulse_ref": da_valc_nr_children_facing_extra_hazard_Reference_Pulse_ref,
            "valc_total_nr_children_facing_extra_hazard_Reference_Pulse": da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse,
            "valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref": da_valc_total_nr_children_facing_extra_hazard_Reference_Pulse_ref,
            "valc_slope_exposure": da_valc_slope_exposure,
            "valc_slope_exposure_ref": da_valc_slope_exposure_ref
        }
    )

    # -------------------------------- Save as Pickles ---------------------------------- #
    with open(data_dir+'{}/source2suffering/reference_pulse/ds_S2S_Reference_Pulse_gmt_{}_{}.pkl'.format(flags['version'],flags['gmt'],flags['rm']), 'wb') as f:
        pk.dump(ds_S2S_Reference_Pulse,f)

    with open(data_dir+'{}/source2suffering/reference_pulse/ds_S2S_Reference_Pulse_Regions_gmt_{}_{}.pkl'.format(flags['version'],flags['gmt'],flags['rm']), 'wb') as f:
        pk.dump(ds_S2S_Reference_Pulse_Regions,f)














#-------------------------------------------------------------------------------------------------------------------------------#
#                                      Source2Suffering Applications for Assessment Reports                                     #
#-------------------------------------------------------------------------------------------------------------------------------#

def assessment_Neptun_Deep():

    print("\n ---------------------------------------------------------")
    print("|          Assessment of the Neptun Deep project          |")
    print("|                for Greenpeace Romania                   |")
    print(" ---------------------------------------------------------")

    # -------------------------------------------------------------------------- #
    # Load the data  computed for the original Expert Report of Wim Thiery       #
    # These are the data that with wich we will compare the results of the       #
    # Source2Suffering framework and its different configurations                #
    # -------------------------------------------------------------------------- #

    # ----------------------------------------------------------- #
    #        valc_nr_extra people facing additionnal hazard       #
    # ----------------------------------------------------------- #

    # Definition of the values per birth cohort

    wt_valc_nr_children_facing_extra_wildfire = [1000]*11
    wt_valc_nr_children_facing_extra_cropfailure = [3000]*7 + [2000]*4
    wt_valc_nr_children_facing_extra_drought = [4000] + [3000]*6 + [2000]*4
    wt_valc_nr_children_facing_extra_flood = [np.nan] * 11
    wt_valc_nr_children_facing_extra_heatwavedarea = [127000,123000,120000,116000,113000,110000,106000,103000,99000,95000,91000]
    wt_valc_nr_children_facing_extra_tropicalcyclone = [1000]*9 + [0]*2

    # Definition of the values for the total birth cohorts of interest

    wt_total_valc_nr_children_facing_extra_wildfire = [11000]
    wt_total_valc_nr_children_facing_extra_cropfailure = [29000]
    wt_total_valc_nr_children_facing_extra_drought = [31000]
    wt_total_valc_nr_children_facing_extra_flood = [np.nan]
    wt_total_valc_nr_children_facing_extra_heatwavedarea = [1203000]
    wt_total_valc_nr_children_facing_extra_tropicalcyclone = [9000]

    # Definition of the values for the total reference birth cohorts of interest

    wt_total_valc_nr_children_facing_extra_wildfire_ref = [0]
    wt_total_valc_nr_children_facing_extra_cropfailure_ref = [0]
    wt_total_valc_nr_children_facing_extra_drought_ref = [0]
    wt_total_valc_nr_children_facing_extra_flood_ref = [np.nan]
    wt_total_valc_nr_children_facing_extra_heatwavedarea_ref = [78000]
    wt_total_valc_nr_children_facing_extra_tropicalcyclone_ref = [0]

    # Creating DataArrays to store the objets  
    
    data = np.array([
        wt_valc_nr_children_facing_extra_wildfire,
        wt_valc_nr_children_facing_extra_cropfailure,
        wt_valc_nr_children_facing_extra_drought,
        wt_valc_nr_children_facing_extra_flood,
        wt_valc_nr_children_facing_extra_heatwavedarea,
        wt_valc_nr_children_facing_extra_tropicalcyclone
    ])

    data = np.flip(np.array([
        wt_valc_nr_children_facing_extra_wildfire,
        wt_valc_nr_children_facing_extra_cropfailure,
        wt_valc_nr_children_facing_extra_drought,
        wt_valc_nr_children_facing_extra_flood,
        wt_valc_nr_children_facing_extra_heatwavedarea,
        wt_valc_nr_children_facing_extra_tropicalcyclone
    ]), axis=1)

    data_total = np.array([
        wt_total_valc_nr_children_facing_extra_wildfire,
        wt_total_valc_nr_children_facing_extra_cropfailure,
        wt_total_valc_nr_children_facing_extra_drought,
        wt_total_valc_nr_children_facing_extra_flood,
        wt_total_valc_nr_children_facing_extra_heatwavedarea,
        wt_total_valc_nr_children_facing_extra_tropicalcyclone    
    ])

    data_total_ref = np.array([
        wt_total_valc_nr_children_facing_extra_wildfire_ref,
        wt_total_valc_nr_children_facing_extra_cropfailure_ref,
        wt_total_valc_nr_children_facing_extra_drought_ref,
        wt_total_valc_nr_children_facing_extra_flood_ref,
        wt_total_valc_nr_children_facing_extra_heatwavedarea_ref,
        wt_total_valc_nr_children_facing_extra_tropicalcyclone_ref    
    ])

    data_total = data_total.squeeze()
    data_total_ref = data_total_ref.squeeze()

    year_start_as = 2010
    year_end_as = 2020

    birth_cohort_int = list(range(year_start_as, year_end_as + 1))  

    hazards = [
        'burntarea', 
        'cropfailedarea', 
        'driedarea', 
        'floodedarea', 
        'heatwavedarea', 
        'tropicalcyclonedarea'
    ]

    da_wt_valc_nr_children_facing_extra_hazard_NeptunDeep = xr.DataArray(
        data,
        dims=["hazard", "birth_year"],
        coords={"hazard": hazards, "birth_year": birth_cohort_int},
        name="wt_valc_nr_children_facing_extra_hazard_NeptunDeep"
    )

    da_wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep = xr.DataArray(
        data_total,
        dims=["hazard"],
        coords={"hazard": hazards},
        name="wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep"
    )

    da_wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep_ref = xr.DataArray(
        data_total_ref,
        dims=["hazard"],
        coords={"hazard": hazards},
        name="wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep_ref"
    )

    # ----------------------------------------------------------- #
    #                    valc_slope_exposure                      #
    # ----------------------------------------------------------- #

    # Importing the .mat objets from WT 
    from scipy.io import loadmat

    slope_exposure_NeptunDeep_heatwave = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_heatwave.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_heatwave = slope_exposure_NeptunDeep_heatwave['slope_exposure_NeptunDeep_heatwave'][:-1]

    slope_exposure_NeptunDeep_heatwave_ref = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_heatwave_ref.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_heatwave_ref = slope_exposure_NeptunDeep_heatwave_ref['slope_exposure_NeptunDeep_heatwave_ref'][:-1]

    slope_exposure_NeptunDeep_cropfailure = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_cropfailure.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_cropfailure = slope_exposure_NeptunDeep_cropfailure['slope_exposure_NeptunDeep_cropfailure'][:-1]

    slope_exposure_NeptunDeep_cropfailure_ref = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_cropfailure_ref.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_cropfailure_ref = slope_exposure_NeptunDeep_cropfailure_ref['slope_exposure_NeptunDeep_cropfailure_ref'][:-1]

    slope_exposure_NeptunDeep_drought = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_drought.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_drought = slope_exposure_NeptunDeep_drought['slope_exposure_NeptunDeep_drought'][:-1]

    slope_exposure_NeptunDeep_drought_ref = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_drought_ref.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_drought_ref = slope_exposure_NeptunDeep_drought_ref['slope_exposure_NeptunDeep_drought_ref'][:-1]

    slope_exposure_NeptunDeep_riverflood = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_riverflood.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_riverflood = slope_exposure_NeptunDeep_riverflood['slope_exposure_NeptunDeep_riverflood'][:-1]

    slope_exposure_NeptunDeep_riverflood_ref = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_riverflood_ref.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_riverflood_ref = slope_exposure_NeptunDeep_riverflood_ref['slope_exposure_NeptunDeep_riverflood_ref'][:-1]

    slope_exposure_NeptunDeep_tropcyclone = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_tropcyclone.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_tropcyclone = slope_exposure_NeptunDeep_tropcyclone['slope_exposure_NeptunDeep_tropcyclone'][:-1]

    slope_exposure_NeptunDeep_tropcyclone_ref = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_tropcyclone_ref.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_tropcyclone_ref = slope_exposure_NeptunDeep_tropcyclone_ref['slope_exposure_NeptunDeep_tropcyclone_ref'][:-1]

    slope_exposure_NeptunDeep_wildfire = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_wildfire.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_wildfire = slope_exposure_NeptunDeep_wildfire['slope_exposure_NeptunDeep_wildfire'][:-1]

    slope_exposure_NeptunDeep_wildfire_ref = loadmat(scripts_dir+'/references/lifetime_exposure_wim/lifetime_exposure_wim_v1/slope_exposure_NeptunDeep_wildfire_ref.mat',squeeze_me=True)
    slope_exposure_NeptunDeep_wildfire_ref = slope_exposure_NeptunDeep_wildfire_ref['slope_exposure_NeptunDeep_wildfire_ref'][:-1]

    # Mapping from extreme types to their corresponding slope exposure arrays
    extreme_to_var = {
        'burntarea': slope_exposure_NeptunDeep_wildfire,
        'cropfailedarea': slope_exposure_NeptunDeep_cropfailure,
        'driedarea': slope_exposure_NeptunDeep_drought,
        'floodedarea': slope_exposure_NeptunDeep_riverflood,
        'heatwavedarea': slope_exposure_NeptunDeep_heatwave,
        'tropicalcyclonedarea': slope_exposure_NeptunDeep_tropcyclone
    }

    # Same mapping for the "_ref" slope exposure arrays (reference scenario)
    extreme_to_var_ref = {
        'burntarea': slope_exposure_NeptunDeep_wildfire_ref,
        'cropfailedarea': slope_exposure_NeptunDeep_cropfailure_ref,
        'driedarea': slope_exposure_NeptunDeep_drought_ref,
        'floodedarea': slope_exposure_NeptunDeep_riverflood_ref,
        'heatwavedarea': slope_exposure_NeptunDeep_heatwave_ref,
        'tropicalcyclonedarea': slope_exposure_NeptunDeep_tropcyclone_ref
    }

    # Build the main DataArray for valc_slope_exposure
    da_valc_slope_exposure = xr.DataArray(
        data=[extreme_to_var[ext] for ext in hazards],
        coords={
            "hazard": hazards,
            "birth_year": birth_cohort_int
        },
        dims=["hazard", "birth_year"],
        name="valc_slope_exposure"
    )

    year_start_as_ref = 1960
    year_end_as_ref = 1970

    birth_cohort_int_ref = list(range(year_start_as_ref, year_end_as_ref + 1))

    # Build the reference DataArray for valc_slope_exposure_ref
    da_valc_slope_exposure_ref = xr.DataArray(
        data=[extreme_to_var_ref[ext] for ext in hazards],
        coords={
            "hazard": hazards,
            "birth_year": birth_cohort_int_ref
        },
        dims=["hazard", "birth_year"],
        name="valc_slope_exposure_ref"
    )

    # ----------------------------------------------------------- #
    #                  ds_WT_NeptunDeep saving                    #
    # ----------------------------------------------------------- #

    # Save as DataSet

    ds_WT_NeptunDeep = xr.Dataset(
    {
        "wt_valc_nr_children_facing_extra_hazard_NeptunDeep": da_wt_valc_nr_children_facing_extra_hazard_NeptunDeep,
        "wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep": da_wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep,
        "wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep_ref" : da_wt_total_valc_nr_children_facing_extra_hazard_NeptunDeep_ref
    }
    )

    # Add both variables to the existing dataset
    ds_WT_NeptunDeep = ds_WT_NeptunDeep.assign(
        valc_slope_exposure=da_valc_slope_exposure,
        valc_slope_exposure_ref=da_valc_slope_exposure_ref
    )

    # Save the DataSet reference as pickles 
    with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_WT_NeptunDeep.pkl'.format(flags['version']), 'wb') as f:
        pk.dump(ds_WT_NeptunDeep,f)

    # -------------------------------------------------------------------------- #
    # Define Total emissions of the fossiel fuel project under study             #
    # -------------------------------------------------------------------------- #
    
    # TOTAL emissions UK oil and gas fields (MtCO2e: add E6 to express as tCO2e) #
    CO2_emissions_NeptunDeep = 207e6

    # -------------------------------------------------------------------------- #
    # Define Transient Climate response to cumulative emission                   #
    # -------------------------------------------------------------------------- #

    # Expert advice (27/03/2024): Je kan 1.65°C per 1000 PgC gebruiken, maar 
    # (mocht dat nuttig zijn) kan je ook transparant gecommuniceerde hogere 
    # percentielen gebruiken. Die geven namelijk een risicoperspectief. 
    # Bijvoorbeeld, 2.2°C per 1000 PgC is nog steeds enkel het 66ste 
    # percentiel. Dat is niet per se extreem. Deze methode is wel vooral 
    # toepasbaar op CO2 and niet CO2-eq. Dus als er veel methaan in die 
    # emissies zit dan zou je dat een beetje moeten aanpassen.
    
    TCRE_init = 0.45 / 1e12 # 1.65°C / 3.7 = 0,45°C per 1000 Gt CO2eq
    
    # -------------------------------------------------------------------------- #
    # GMT produce by Neptun Deep using the TCRE                                  #
    # -------------------------------------------------------------------------- #

    print("\nCO2 emissions of the Neptun Deep Project = {} MtCO2eq\n".format(CO2_emissions_NeptunDeep/10**6))

    print("Value used for TCRE = 0.45 °C per 1000 Gt CO2eq\n")

    print("GMT produce by Neptun Deep = {} °C".format(TCRE_init*CO2_emissions_NeptunDeep))

    # -------------------------------------------------------------------------- #
    # Question 1: call function to compute the number of children in the world   #
    # born between 2010 and 2020 facing one extra heatwave due to the total      # 
    # emissions of the three fields of Neptun Deep                               #
    # -------------------------------------------------------------------------- #

    print("")
    print(" -----------------------------------------------------------------------------------------------------------------")
    print("| Q.(1) - Number of children in the world facing an additional heatwave due to the total emissions of Neptun Deep |")
    print(" -----------------------------------------------------------------------------------------------------------------")
    print("")

    if flags['old_demo']:

        # Load the corresponding exposure dataset with old demography
        with open(data_dir+'{}/old_demography/{}/ds_le_perregion_gmt_{}_{}.pkl'.format('pickles_sandbox','heatwavedarea',flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_perregion = pk.load(f)

        # Load the absolute cohort sizes at the regional level with the old demography
        with open(data_dir + '{}/old_demography/country/da_valp_cohort_size_abs.pkl'.format('pickles_sandbox'), 'rb') as f:
            da_valp_cohort_size_abs = pk.load(f)

    else: 

        # Load the corresponding exposure dataset with new demography
        with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],'heatwavedarea',flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_perregion = pk.load(f)

        # Load the absolute cohort sizes at the regional level with the new demography
        with open(data_dir + '{}/country/da_valp_cohort_size_abs.pkl'.format(flags['version']), 'rb') as f:
            da_valp_cohort_size_abs = pk.load(f)

    # Computation for the 2010 to 2020 birth cohorts #

    year_start_as = 2010
    year_end_as = 2020

    valc_nr_children_facing_extra_heatwave_NeptunDeep, S2S_slope_exposure = emissions2npeople(
        CO2_emissions = CO2_emissions_NeptunDeep,
        TCRE = TCRE_init,
        ds_le = ds_le_perregion,
        region_ind = 11,
        year_start = year_start_as,
        year_end = year_end_as,
        df_GMT_strj = df_GMT_strj, 
        da_valp_cohort_size_abs = da_valp_cohort_size_abs,
        rounding = 2)  
    
    # Generate list of birth years for iteration
    years_loop = list(range(year_end_as, year_start_as - 1, -1))
    nbirthyears = len(years_loop)

    for i in range(nbirthyears):
    
        print("For birth year {} = {} children".format(years_loop[i], int(valc_nr_children_facing_extra_heatwave_NeptunDeep[i])))
    
    print("For total birth of the {}-{} period = {} children \n".format(year_start_as,year_end_as, int(valc_nr_children_facing_extra_heatwave_NeptunDeep[-1])))

    # Computation for the reference 1960 - 1970 birth cohorts # 

    year_start_as_ref = 1960
    year_end_as_ref = 1970

    # Generate list of birth years for iteration
    #years_loop = list(range(year_end_as_ref, year_start_as_ref - 1, -1))

    valc_nr_children_facing_extra_heatwave_NeptunDeep_ref, S2S_slope_exposure_ref = emissions2npeople(
        CO2_emissions = CO2_emissions_NeptunDeep,
        TCRE = TCRE_init,
        ds_le = ds_le_perregion,
        region_ind = 11,
        year_start = year_start_as_ref,
        year_end = year_end_as_ref,
        df_GMT_strj = df_GMT_strj, 
        da_valp_cohort_size_abs = da_valp_cohort_size_abs,
        rounding = 2)

    # Only keep the value for the total of the birth cohort between 1960 and 1970

    valc_nr_children_facing_extra_heatwave_NeptunDeep_ref = valc_nr_children_facing_extra_heatwave_NeptunDeep_ref[-1]

    # Compute the relative exposure to an additionnal heatwaves between the 2010-2020 and the 1960-1970 birth cohort
    
    if valc_nr_children_facing_extra_heatwave_NeptunDeep_ref != 0:

        rel_percentage = round((valc_nr_children_facing_extra_heatwave_NeptunDeep[-1] / valc_nr_children_facing_extra_heatwave_NeptunDeep_ref) * 100)
    
    else:

        rel_percentage = np.nan

    valc_nr_children_facing_extra_heatwave_NeptunDeep_ref = [
        valc_nr_children_facing_extra_heatwave_NeptunDeep_ref,
        rel_percentage
    ]

    print("Number of people born between 1960 and 1970 that will be exposed to an additionnal heatwave due to the total emissions of Neptun Deep = {} people".format(int(valc_nr_children_facing_extra_heatwave_NeptunDeep_ref[0])))
    print("Relative difference between the 2010-2020 and the 1960-1970 birth cohorts = {} %".format(valc_nr_children_facing_extra_heatwave_NeptunDeep_ref[-1]))

    # # -------------------------------------------------------------------------- #
    # # Question 2: call function to compute the number of people facing one       #
    # # extra ... due to emissions                                                 #
    # # -------------------------------------------------------------------------- #

    print("")
    print(" ---------------------------------------------------------------------------------------------------------------")
    print("| Q.(2) - Number of children in the world facing an additional hazard due to the total emissions of Neptun Deep |")
    print(" ---------------------------------------------------------------------------------------------------------------")
    print("")

    year_start_as = 2010
    year_end_as = 2020
    year_start_as_ref = 1960
    year_end_as_ref = 1970
    
    all_extremes = [
        'burntarea', 
        'cropfailedarea', 
        'driedarea', 
        'floodedarea', 
        'heatwavedarea', 
        'tropicalcyclonedarea'
    ]

    hazards_name = [
        'Wild Fires', 
        'Crop Failures', 
        'Droughts', 
        'River Floods', 
        'Heatwaves', 
        'Tropical Cyclones'
    ]

    # Dictionary to store values for each extreme hazard type
    d_valc = {}
    d_valc_ref = {}
    d_valc_total = {}
    d_valc_total_ref = {}
    d_valc_slope_exposure = {}
    d_valc_slope_exposure_ref = {}

    n = 0
    # Loop over all extreme hazard types
    for extr in all_extremes:

        extr_name = hazards_name[n]
        print("---------- Hazard = {} ----------\n".format(extr_name))
        n+=1

        if flags['old_demo']:

            # Load the corresponding exposure dataset with old demography
            with open(data_dir+'{}/old_demography/{}/ds_le_perregion_gmt_{}_{}.pkl'.format('pickles_sandbox',extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_perregion = pk.load(f)

            # Load the absolute cohort sizes at the regional level with the old demography
            with open(data_dir + '{}/old_demography/country/da_valp_cohort_size_abs.pkl'.format('pickles_sandbox'), 'rb') as f:
                da_valp_cohort_size_abs = pk.load(f)

        else: 

            # Load the corresponding exposure dataset with new demography
            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_perregion = pk.load(f)

            # Load the absolute cohort sizes at the regional level with the new demography
            with open(data_dir + '{}/country/da_valp_cohort_size_abs.pkl'.format(flags['version']), 'rb') as f:
                da_valp_cohort_size_abs = pk.load(f)

        # Compute the number of children exposed to the extra hazard under the NeptunDeep scenario
        valc_nr_children_facing_extra_hazard_NeptunDeep, S2S_slope_exposure = emissions2npeople(
            CO2_emissions=CO2_emissions_NeptunDeep,
            TCRE=TCRE_init,
            ds_le=ds_le_perregion,
            region_ind=11,
            year_start=year_start_as,
            year_end=year_end_as,
            df_GMT_strj=df_GMT_strj,
            da_valp_cohort_size_abs=da_valp_cohort_size_abs,
            rounding=2
        )

        # Compute the number of children exposed to the extra hazard under the NeptunDeep scenario for the reference birth cohorts
        valc_nr_children_facing_extra_hazard_NeptunDeep_ref, S2S_slope_exposure_ref = emissions2npeople(
            CO2_emissions=CO2_emissions_NeptunDeep,
            TCRE=TCRE_init,
            ds_le=ds_le_perregion,
            region_ind=11,
            year_start=year_start_as_ref,
            year_end=year_end_as_ref,
            df_GMT_strj=df_GMT_strj,
            da_valp_cohort_size_abs=da_valp_cohort_size_abs,
            rounding=2
        )

        nr_total = valc_nr_children_facing_extra_hazard_NeptunDeep[-1]
        nr_total_ref = valc_nr_children_facing_extra_hazard_NeptunDeep_ref[-1]

        # Store the value of the total number of people exposed
        d_valc_total[extr] = nr_total
        d_valc_total_ref[extr] = nr_total_ref

        # Store the value of slope exposure for each hazard

        d_valc_slope_exposure[extr] = S2S_slope_exposure
        d_valc_slope_exposure_ref[extr] = S2S_slope_exposure_ref

        # Exclude the last value along the only dimension (since the array is 1D)
        da_trimmed = valc_nr_children_facing_extra_hazard_NeptunDeep[:-1]
        da_ref_trimmed = valc_nr_children_facing_extra_hazard_NeptunDeep_ref[:-1]

        # Reverse the values
        da_reversed = da_trimmed[::-1]
        da_ref_reversed = da_ref_trimmed[::-1]

        # Store the result in the dictionary under the key 'extr'
        d_valc[extr] = da_reversed
        d_valc_ref[extr] = da_ref_reversed

        # Generate list of birth years in descending order
        years_loop = list(range(year_end_as, year_start_as - 1, -1))
        nbirthyears = len(years_loop)
        years_loop_ref = list(range(year_end_as_ref, year_start_as_ref - 1, -1))
        nbirthyears_ref = len(years_loop_ref)

        # Print number of exposed children per birth year
        for i in range(nbirthyears):
            print("For birth year {} = {} children".format(
                years_loop[i], int(valc_nr_children_facing_extra_hazard_NeptunDeep[i])))

        # Print total number of exposed children across all birth years
        print("For total birth of the {}–{} period = {} children \n".format(
            year_start_as, year_end_as, int(valc_nr_children_facing_extra_hazard_NeptunDeep[-1]))
        )
        # Print total number of exposed children across all birth years for the reference period
        print("For total birth of the {}–{} period = {} children \n".format(
            year_start_as_ref, year_end_as_ref, int(valc_nr_children_facing_extra_hazard_NeptunDeep_ref[-1]))
        )

    # ----------------------- Creating DataArray and DataSet ---------------------------- #

    # Create a DataArray containing all values for each hazard and birth year
    da_valc_nr_children_facing_extra_hazard_NeptunDeep = xr.DataArray(
        data=[d_valc[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as, year_end_as + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_nr_children_facing_extra_hazard_NeptunDeep"
    )

    da_valc_nr_children_facing_extra_hazard_NeptunDeep_ref = xr.DataArray(
        data=[d_valc_ref[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as_ref, year_end_as_ref + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_nr_children_facing_extra_hazard_NeptunDeep_ref"
    )

    # Create a DataArray containing the values of the total of people exposed for each hazard
    da_valc_total_nr_children_facing_extra_hazard_NeptunDeep = xr.DataArray(
        data=[d_valc_total[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
        },
        dims=["hazard"],
        name="valc_total_nr_children_facing_extra_hazard_NeptunDeep"
    )

    da_valc_total_nr_children_facing_extra_hazard_NeptunDeep_ref = xr.DataArray(
        data=[d_valc_total_ref[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
        },
        dims=["hazard"],
        name="valc_total_nr_children_facing_extra_hazard_NeptunDeep_ref"
    )

    # Create a DataArray containing the values of the slope exposure for each hazard

    da_valc_slope_exposure = xr.DataArray(
        data=[d_valc_slope_exposure[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as, year_end_as + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_slope_exposure"
    )

    da_valc_slope_exposure_ref = xr.DataArray(
        data=[d_valc_slope_exposure_ref[extr] for extr in all_extremes],
        coords={
            "hazard": all_extremes,
            "birth_year": list(range(year_start_as_ref, year_end_as_ref + 1))  # in ascending order
        },
        dims=["hazard", "birth_year"],
        name="valc_slope_exposure_ref"
    )

    ds_S2S_NeptunDeep = xr.Dataset(
    {
        "valc_nr_children_facing_extra_hazard_NeptunDeep": da_valc_nr_children_facing_extra_hazard_NeptunDeep,
        "valc_total_nr_children_facing_extra_hazard_NeptunDeep": da_valc_total_nr_children_facing_extra_hazard_NeptunDeep,
        "valc_slope_exposure": da_valc_slope_exposure,
        "valc_nr_children_facing_extra_hazard_NeptunDeep_ref": da_valc_nr_children_facing_extra_hazard_NeptunDeep_ref,
        "valc_total_nr_children_facing_extra_hazard_NeptunDeep_ref": da_valc_total_nr_children_facing_extra_hazard_NeptunDeep_ref,
        "valc_slope_exposure_ref": da_valc_slope_exposure_ref

    }
    )

    
    # -------------------------------- Save as Pickles ---------------------------------- #

    if flags['old_demo']:

        # dump pickle of ds_valc_nr_children_facing_extra_hazard_NeptunDeep with the old demography
        with open(data_dir+'{}/old_demography/assessment/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format('pickles_sandbox',flags['gmt'],flags['rm']), 'wb') as f:
            pk.dump(ds_S2S_NeptunDeep,f)

    else: 

        # dump pickle of ds_valc_nr_children_facing_extra_hazard_NeptunDeep with the new demography
        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],flags['gmt'],flags['rm']), 'wb') as f:
            pk.dump(ds_S2S_NeptunDeep,f)
    
    # -------------------------------------------------------------------------- #
    # Question 3: compute the number of heat-related deaths between              #
    # today and 2100 due to the total emissions of the three fields              #
    # -------------------------------------------------------------------------- #

    valc_mortality_NeptunDeep = mortality_cost_carbon(CO2_emissions=CO2_emissions_NeptunDeep)
    
    print("")
    print(" ------------------------------------------------------------------------------------------------------------------------")
    print("| Q.(3) - Number of heat-related deaths that are expected worldwide until 2100 due to the total emissions of Neptun Deep |")
    print(" ------------------------------------------------------------------------------------------------------------------------")
    print("")


    print("Mortality Cost of Carbon = {} heat-related deaths until 2100".format(valc_mortality_NeptunDeep))


def assessment_Nordic_Barents_Sea():

    print("\n ---------------------------------------------------------")
    print("|          Assessment of the Barents Sea project          |")
    print("|                for Greenpeace Nordic                    |")
    print(" ---------------------------------------------------------")

    # Still must be finish #