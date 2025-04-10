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
from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

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

        

        # if flags['global_avg_emergence']:
        #     print("---------------------------------------------------------")
        #     print("Report 1 : Estimates of land area for 1960 and 2020 emergence of multiple hazards")
        #     # estimates of land area and (potential) pf for 1960 and 2020 emergencve of multiple hazards
        #     multi_hazard_emergence(
        #         grid_area,
        #         da_emergence_mean,
        #         da_gs_popdenom,
        #     )
        # else:
        #     print("---------------------------------------------------------")
        #     print("Report 1 not produced because flags['global_avg_emergence']=0")
        
        # # get birth year cohort sizes at grid scale
        # print("---------------------------------------------------------")
        # print("Report 2 : Gridscale cohort sizes per birth year.")
        # gridscale_cohort_sizes(
        #     flags,
        #     da_population,
        #     gridscale_countries,   
        # )    
        # print("Save gridscale_cohort_global.pkl in data/{version}/country")
        
        # # per hazard, locations where exposure occurs across whole ensemble
        # print("---------------------------------------------------------")
        # print("Report 3 : Grid exposure locations for all simulations.")
        # exposure_locs(
        #     flags,
        #     grid_area,
        # )
        # print("Save exposure_occurrence_{extr}.pkl in data/{version}/{extr}")
        
        # # per run for 1.5, 2.5, 2.7 and 3.5, collect maps of emergence locations to be used in geographically constrained pf estimates
        # print("---------------------------------------------------------")
        # print("Report 4: Collect maps of emergence locations")
        # emergence_locs_perrun(
        #     flags,
        #     grid_area,
        #     gridscale_countries,
        #     countries_mask,
        #     countries_regions,
        # )    
        # print("Save emergence_locs_perrun_{extr}_{step}.pkl in /data/{version}/{extr}")

        # # compute geographically constrained pf
        # print("---------------------------------------------------------")
        # print("Report 5: Population fraction estimates per run for selected GMTs when constraining denominator by geography")
        # pf_geoconstrained(
        #     flags,
        #     countries_mask,
        # )
        # print("Save pf_geoconstrained_{extr}.pkl in data/{version}/{extr}")
        
        # # print geographically constrained pf vs regular pf
        # print("---------------------------------------------------------")
        # print("Report 6: Print geographically constrained pf Vs regular pf")
        # print_pf_geoconstrained(
        #     flags,
        #     da_gs_popdenom,
        # )    

        # # checking for signifiance of change in means between 1960 and 2020 pf per event and for a GMT level
        # print("---------------------------------------------------------")
        # print("Report 7: Check for signifiance of change in means between 1960 and 2020 pf per event and for a GMT level")
        # print("Not execute due to bugs in the original script")
        # # paired_ttest(
        # #     flags,
        # #     da_gs_popdenom,
        # # )
        
        # # print latex table on ensemble members per hazard
        # print("---------------------------------------------------------")
        # print("Report 8: Print LaTeX table on ensemble members per hazard")
        # if latex_output == True:
        #     print_latex_table_ensemble_sizes(
        #         flags,
        #         df_GMT_strj,
        #     )
        # else:
        #     print('No LaTeX output in the configuration (see reporting.py)')   

        # # children (i.e. those born between 2003-2020) living unprec exposure between 1.5 and 2.7 degrees warming (for numbers in conclusion of paper)
        # print("---------------------------------------------------------")    
        # print("Report 9: Number of children (born between 2003-2020) living ULE between 1.5 and 2.7°C warming")
        # print_millions_excess(
        #     flags,
        #     df_GMT_strj,
        # )     

        # # print pf info 
        # print("---------------------------------------------------------") 
        # print("Report 10: Ratio of pfs")
        # print_pf_ratios_and_abstract_numbers(
        #     flags,
        #     df_GMT_strj,
        #     da_gs_popdenom,
        # )    
        
        # # get number of million people unprecedented: (will change this stuff to run for all extremes and birth years for table in paper)
        # print("---------------------------------------------------------")
        # print("Report 11: Number of million people unprecedented")
        # print("Not execute due to bugs in the original script")
        # # print_absolute_unprecedented(
        # #     ds_pf_gs,
        # # )
    
        # get cities that work for f1 concept plot
        print("---------------------------------------------------------")
        print("Report 12: Get pickle of cities that are valid for f1 concept plot") 
        find_valid_cities(
            df_countries,
            da_cohort_size,
            countries_mask,
            countries_regions,
            d_isimip_meta,
            flags,
        )
        
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
    
    print("-------------------------------------------------------------------------")
    print("Start source2suffering_report framework from the Source2Suffering Project")
    print("-------------------------------------------------------------------------")

    sys.path.append(os.path.abspath(scripts_dir+"/output/papers/source2suffering"))
    from source2suffering_report import *

    # Configuration of the reports #

    validation_backward_comp = 0        # 0: do not produce reports for backward compatibility of 2025 Python Lifetime Exposure Framework with 2021 Matlab Thiery Lifetime Exposure Framework
                                        # 1: produce reports for backward compatibility of 2025 Python Lifetime Exposure Framework with 2021 Matlab Thiery Lifetime Exposure Framework

    if validation_backward_comp:

        valRomania = 1                  # 0: do not produce reports for Neptun Deep project
                                        # 1: produce reports for Neptun Deep project
        valNorLic = 0                   # 0: do not produce reports for testimony in Norwegian lawsuit
                                        # 1: produce reports for testimony in Norwegian lawsuit
        valNorAppea = 0                 # 0: do not produce reports for Norwegian lawsuit - written report June 2024 for appeal
                                        # 1: produce reports for Norwegian lawsuit - written report June 2024 for appeal
        valNorECtHR = 0                 # 0: do not produce reports for Norwegian lawsuit - written report ECtHR June 2024 for appeal
                                        # 1: produce reports for Norwegian lawsuit - written report ECtHR June 2024 for appeal

    if valRomania:

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
        
        TCRE = 0.45 / 1e12 # 1.65°C / 3.7 = 0,45°C per 1000 Gt CO2eq

        # -------------------------------------------------------------------------- #
        # Define Mortality cost of Carbon                                            #
        # -------------------------------------------------------------------------- #

        # 4,434 metric tons of carbon dioxide in 2020 cfr. (Bressler, 2021) 
        # causes one excess death globally in expectation between 2020-2100 
        mortality_cost_carbon = 4434 
        
        # -------------------------------------------------------------------------- #
        # Question 1: call function to compute the number of people facing one       #
        # extra heatwave due to the total emissions of the three fields              #
        # -------------------------------------------------------------------------- #

        valc_nr_children_facing_extra_heatwave_NeptunDeep = emissions2npeople(
            CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, 
            GMT_BE, valp_cohort_size_abs, 1, 5)  # 5 is for heatwave

        # Include a small add-on table for the reference years 1960-1970

        valc_nr_children_facing_extra_heatwave_NeptunDeep_ref = emissions2npeople(
            CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, 
            GMT_BE, valp_cohort_size_abs, 1, 5)  # 5 is for heatwaves

        valc_nr_children_facing_extra_heatwave_NeptunDeep_ref = valc_nr_children_facing_extra_heatwave_NeptunDeep_ref[-1]

        valc_nr_children_facing_extra_heatwave_NeptunDeep_ref = [
            valc_nr_children_facing_extra_heatwave_NeptunDeep_ref,
            round(valc_nr_children_facing_extra_heatwave_NeptunDeep[-1] / valc_nr_children_facing_extra_heatwave_NeptunDeep_ref * 100)
        ]

        # -------------------------------------------------------------------------- #
        # Question 2: call function to compute the number of people facing one       #
        # extra ... due to emissions                                                 #
        # -------------------------------------------------------------------------- #

        impacts = {1: "wildfire", 2: "crop failure", 3: "drought", 4: "river floods", 5: "heatwaves", 6: "tropical cyclones"}
        valc_impacts_NeptunDeep = {}
        valc_impacts_NeptunDeep_ref = {}

        for key, impact in impacts.items():
            valc_impacts_NeptunDeep[impact] = mf_emissions2npeople(
                CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, 
                GMT_BE, valp_cohort_size_abs, 1, key)
    
            # Include a small add-on table for the reference years 1960-1970

            valc_impacts_NeptunDeep_ref[impact] = mf_emissions2npeople(
                CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, 
                GMT_BE, valp_cohort_size_abs, 1, key)
            valc_impacts_NeptunDeep_ref[impact] = valc_impacts_NeptunDeep_ref[impact][-1]
            valc_impacts_NeptunDeep_ref[impact] = [
                valc_impacts_NeptunDeep_ref[impact],
                round(valc_impacts_NeptunDeep[impact][-1] / valc_impacts_NeptunDeep_ref[impact] * 100)
            ]

        # -------------------------------------------------------------------------- #
        # Question 3: compute the number of heat-related deaths between              #
        # today and 2100 due to the total emissions of the three fields              #
        # -------------------------------------------------------------------------- #

        valc_mortality_NeptunDeep = (CO2_emissions_NeptunDeep / mortality_cost_carbon // 1000) * 1000



#%%-------------------------------------------------------------------------------------------#
# Framework to report not configured                                                          #
#---------------------------------------------------------------------------------------------#

if not(Thiery_2021 or Grant_2025 or Laridon_2025 or Source2Suffering):
    print("No pre-defined Report Framework outside Thiery et al.(2021), Grant et al.(2025), Laridon et al.(2025) or Source2Suffering Project.")
