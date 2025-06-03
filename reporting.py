# ----------------------------------------------------------------------------------------#
# Subscript to execute the computations and outputs of the different reports              #
# needed in the associate papers                                                          #
# ----------------------------------------------------------------------------------------#
#%%-------------------------------------------------------------- #
# Libraries                                                       #
# --------------------------------------------------------------- #

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

#%%-----------------------------------------------------------------------------------#
# Framework to compute and print all the reports associated to Thiery et al.(2021)    #
# This framework has not been reproduced for the moment.                              #
#-------------------------------------------------------------------------------------#

if Thiery_2021:

    print("-----------------------------------------------------------")
    print("Start thiery_2021_report framework from Thiery et al.(2021)")
    print("-----------------------------------------------------------")

    sys.path.append(os.path.abspath(scripts_dir+"/output/papers/thiery_2021"))
    from thiery_2021_report import *

    pass

#%%-----------------------------------------------------------------------------------#
# Framework to compute and print all the reports associated to Grant et al.(2025)     #
#-------------------------------------------------------------------------------------#

if Grant_2025:
 
    paper = True          # produce the reports needed for the paper of Grant et al.(2025)
    latex_output = False  # decide to print the LaTeX code for tables

    if paper : 

        print("---------------------------------------------------------")
        print("Start grant_2025_report framework from Grant et al.(2025)")
        print("---------------------------------------------------------")

        sys.path.append(os.path.abspath(scripts_dir+"/output/papers/grant_2025"))
        from grant_2025_report import *

        if flags['global_avg_emergence']:
            print("---------------------------------------------------------")
            print("Report 1 : Estimates of land area for 1960 and 2020 emergence of multiple hazards")
            # estimates of land area and (potential) pf for 1960 and 2020 emergencve of multiple hazards
            multi_hazard_emergence(
                grid_area,
                da_emergence_mean,
                da_gs_popdenom,
            )
        else:
            print("---------------------------------------------------------")
            print("Report 1 not produced because flags['global_avg_emergence']=0")
        
        # get birth year cohort sizes at grid scale
        print("---------------------------------------------------------")
        print("Report 2 : Gridscale cohort sizes per birth year.")
        gridscale_cohort_sizes(
            flags,
            da_population,
            gridscale_countries,   
        )    
        print("Save gridscale_cohort_global.pkl in data/{version}/country")
        
        # per hazard, locations where exposure occurs across whole ensemble
        print("---------------------------------------------------------")
        print("Report 3 : Grid exposure locations for all simulations.")
        exposure_locs(
            flags,
            grid_area,
        )
        print("Save exposure_occurrence_{extr}.pkl in data/{version}/{extr}")
        
        # per run for 1.5, 2.5, 2.7 and 3.5, collect maps of emergence locations to be used in geographically constrained pf estimates
        print("---------------------------------------------------------")
        print("Report 4: Collect maps of emergence locations")
        emergence_locs_perrun(
            flags,
            grid_area,
            gridscale_countries,
            countries_mask,
            countries_regions,
        )    
        print("Save emergence_locs_perrun_{extr}_{step}.pkl in /data/{version}/{extr}")

        # compute geographically constrained pf
        print("---------------------------------------------------------")
        print("Report 5: Population fraction estimates per run for selected GMTs when constraining denominator by geography")
        pf_geoconstrained(
            flags,
            countries_mask,
        )
        print("Save pf_geoconstrained_{extr}.pkl in data/{version}/{extr}")
        
        # print geographically constrained pf vs regular pf
        print("---------------------------------------------------------")
        print("Report 6: Print geographically constrained pf Vs regular pf")
        print_pf_geoconstrained(
            flags,
            da_gs_popdenom,
        )    

        # checking for signifiance of change in means between 1960 and 2020 pf per event and for a GMT level
        print("---------------------------------------------------------")
        print("Report 7: Check for signifiance of change in means between 1960 and 2020 pf per event and for a GMT level")
        print("Not execute due to bugs in the original script")
        # paired_ttest(
        #     flags,
        #     da_gs_popdenom,
        # )
        
        # print latex table on ensemble members per hazard
        print("---------------------------------------------------------")
        print("Report 8: Print LaTeX table on ensemble members per hazard")
        if latex_output == True:
            print_latex_table_ensemble_sizes(
                flags,
                df_GMT_strj,
            )
        else:
            print('No LaTeX output in the configuration (see reporting.py)')   

        # children (i.e. those born between 2003-2020) living unprec exposure between 1.5 and 2.7 degrees warming (for numbers in conclusion of paper)
        print("---------------------------------------------------------")    
        print("Report 9: Number of children (born between 2003-2020) living ULE between 1.5 and 2.7°C warming")
        print_millions_excess(
            flags,
            df_GMT_strj,
        )     

        # print pf info 
        print("---------------------------------------------------------") 
        print("Report 10: Ratio of pfs")
        print_pf_ratios_and_abstract_numbers(
            flags,
            df_GMT_strj,
            da_gs_popdenom,
        )    
        
        # get number of million people unprecedented: (will change this stuff to run for all extremes and birth years for table in paper)
        print("---------------------------------------------------------")
        print("Report 11: Number of million people unprecedented")
        print("Not execute due to bugs in the original script")
        # print_absolute_unprecedented(
        #     ds_pf_gs,
        # )
    
        # get cities that work for f1 concept plot
        print("---------------------------------------------------------")
        print("Report 12: Get pickle of cities that are valid for f1 concept plot") 
        print("Not execute due to bugs in the original script")
        # find_valid_cities(
        #     df_countries,
        #     da_cohort_size,
        #     countries_mask,
        #     countries_regions,
        #     d_isimip_meta,
        #     flags,
        # )
        
        # latex tables of CF per extr and GMT pathway
        print("---------------------------------------------------------")
        print("Report 13: LaTeX tables of CF per extreme and GMT pathway")
        if latex_output==True:
            print_latex_table_unprecedented(
                flags,
                da_gs_popdenom,
            )
        else:
            print('No LaTeX output in the configuration (see reporting.py)')     
        
        # large latex tables of CF per extr, country and pathway
        print("---------------------------------------------------------")
        print("Report 14: Large LaTeX tables on CF data per country, birth year and 1.5, 2.5 and 3.5 degree scenario")
        print("Not execute due to bugs in the original script")
        if latex_output==True:
            pass
            # print_latex_table_unprecedented_sideways(
            #     flags,
            #     da_gs_popdenom,
            # )
        else:
            print('No LaTeX output in the configuration (see reporting.py)')     
        
        # data for box plots of heatwaves (f2)
        print("---------------------------------------------------------")
        print("Report 15: Data for boxplots of heatwaves in f2")
        print_f2_info(
            ds_pf_gs,
            flags,
            df_GMT_strj,
            da_gs_popdenom,
            gdf_country_borders,
        )
        
        # data for f3
        print("---------------------------------------------------------")
        print("Report 16: Data for f3")
        print_f3_info(
            flags,
            da_gs_popdenom
        )
        
        # data for pyramid stuff (f4)
        print("---------------------------------------------------------")
        print("Report 17: Data for f4")
        print("Not execute because data are missing")
        # print_pyramid_info(
        #     flags,
        # )

#%%-----------------------------------------------------------------------------------#
# Framework to compute and print all the reports associated to the Assessment reports # 
#-------------------------------------------------------------------------------------#

# Assessment based on Grant et al.(2025) version
if Grant_2025:  

    # Configuration of which reports are produced # 

    Norway_BiCC2 = 0      # produce the reports needed for the Norway BiCC2 assessment
    latex_output = False  # decide to print the LaTeX code for tables

    if Norway_BiCC2:
        
            sys.path.append(os.path.abspath(scripts_dir+"/output/assessment"))
            from report_assessment import *

            # latex tables of CF per extr and GMT pathway
            print("-------------------------------------------------------------------")
            print("Report 1: LaTeX tables of CF per extreme and GMT pathway for Norway")
            if latex_output==True:
                print_latex_table_unprecedented(
                    flags,
                    da_gs_popdenom,
                )
            else:
                print('No LaTeX output in the configuration (see reporting.py)')

#%%-----------------------------------------------------------------------------------#
# Framework to compute and print all the reports associated to Laridon et al.(2025)   # 
#-------------------------------------------------------------------------------------#

if Laridon_2025:

    print("--------------------------------------------------------------")
    print("Start laridon_2025_report framework from  Laridon et al.(2025)")
    print("--------------------------------------------------------------")
    
    sys.path.append(os.path.abspath(scripts_dir+"/output/papers/laridon_2025"))
    from laridon_2025_report import *
    
    print('Framework to report for Laridon et al.(2025) not yet configured')

#%%-------------------------------------------------------------------------------------------#
# Framework to compute and print all the reports associated to the Source2Suffering Project   # 
#---------------------------------------------------------------------------------------------#

if Source2Suffering:
    
    # sys.path.append(os.path.abspath(scripts_dir+"/output/papers/source2suffering"))
    # from source2suffering_report import *

    from source2suffering import *

    # Configuration of the reports #

    validation_backward_comp = 1        # 0: do not produce reports for backward compatibility of 2025 Python Lifetime Exposure Framework with 2021 Matlab Thiery Lifetime Exposure Framework
                                        # 1: produce reports for backward compatibility of 2025 Python Lifetime Exposure Framework with 2021 Matlab Thiery Lifetime Exposure Framework

    if validation_backward_comp:

        Greenpeace_Romania_Neptun_Deep = 1      # 0: do not produce reports of Neptun Deep project for Greenpeace Romania
                                                # 1: produce reports of Neptun Deep project for Greenpeace Romania
        Greenpeace_Nordic_Barents_Sea = 0       # 0: do not produce reports of Barents Sea project for Greenpeace Nordic
                                                # 1: produce reports of Barents Sea project for Greenpeace Nordic


    if Greenpeace_Romania_Neptun_Deep:

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

        # Build the reference DataArray for valc_slope_exposure_ref
        da_valc_slope_exposure_ref = xr.DataArray(
            data=[extreme_to_var_ref[ext] for ext in hazards],
            coords={
                "hazard": hazards,
                "birth_year": birth_cohort_int
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

        # Load the corresponding exposure dataset
        with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],'heatwavedarea',flags['gmt'],flags['rm']), 'rb') as f:
            ds_le_perregion = pk.load(f)

        # Load the absolute cohort sizes at the regional level
            with open(data_dir + '{}/country/da_valp_cohort_size_abs.pkl'.format(flags['version'])) as f:
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

            # Load the corresponding exposure dataset
            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_perregion = pk.load(f)

            # Load the absolute cohort sizes at the regional level
            with open(data_dir + '{}/country/da_valp_cohort_size_abs.pkl'.format(flags['version'])) as f:
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

        # dump pickle of ds_valc_nr_children_facing_extra_hazard_NeptunDeep
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

    if Greenpeace_Nordic_Barents_Sea:

        print("\n ---------------------------------------------------------")
        print("|          Assessment of the Barents Sea project          |")
        print("|                for Greenpeace Nordic                    |")
        print(" ---------------------------------------------------------")

#%%-------------------------------------------------------------------------------------------#
# Framework to report not configured                                                          #
#---------------------------------------------------------------------------------------------#

if not(Thiery_2021 or Grant_2025 or Laridon_2025 or Source2Suffering):
    print("No pre-defined Report Framework outside Thiery et al.(2021), Grant et al.(2025), Laridon et al.(2025) or Source2Suffering Project.")
