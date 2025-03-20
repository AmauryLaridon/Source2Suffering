# ----------------------------------------------------------------------------------------
# Script to compute and print the different reports needed in the associate papers
# ----------------------------------------------------------------------------------------
#%%  ----------------------------------------------------------------
# Libraries  
# ----------------------------------------------------------------

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
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins = init()


#%%------------------------------------------------------------------------------------
# Framework to compute and print all the reports associated to Thiery et al.(2021) 
#--------------------------------------------------------------------------------------

if Thiery_2021 == True:
    pass

#%%------------------------------------------------------------------------------------
# Framework to compute and print all the reports associated to Grant et al.(2025) 
#--------------------------------------------------------------------------------------

if Grant_2025 == True:

    print("---------------------------------------------------------")
    print("Start grant_2025_report framework from Grant et al.(2025)")
    print("---------------------------------------------------------")

    sys.path.append(os.path.abspath(scripts_dir+"/output/papers/grant_2025"))
    from grant_2025_report import *

    latex_output = False #decide to print the LaTeX code for tables

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
    if flags['version']=='pickles_v2':
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
    else:
        print("---------------------------------------------------------") 
        print("Report 12: Not execute because needs data from {version} = 'pickles_v2'")
    
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


#%%------------------------------------------------------------------------------------
# Framework to compute and print all the reports associated to Laridon et al.(2025) 
#--------------------------------------------------------------------------------------

if Laridon_2025 == True:
    pass


#%%------------------------------------------------------------------------------------
# Framework to report not configured 
#--------------------------------------------------------------------------------------

if not(Thiery_2021 or Grant_2025 or Laridon_2025):
    print("No pre-defined Report Framework")
