# -----------------------------------------------------------------------------------------------
# Scripts with plots functions and produce the plots based on the configuration of the Framework
# -----------------------------------------------------------------------------------------------

#%%  ----------------------------------------------------------------
# Libraries 
# ----------------------------------------------------------------

import os
import sys
import requests
from zipfile import ZipFile
import io
import xarray as xr
import pickle as pk
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as pe
import mapclassify as mc
from copy import deepcopy as cp
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import regionmask as rm
import geopandas as gpd
from scipy import interpolate
from scipy import stats as sts
import cartopy.crs as ccrs
import seaborn as sns
import cartopy as cr
import cartopy.feature as feature
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

#%%------------------------------------------------------------------------------------
# Framework to plots all figures associated to Thiery et al.(2021)
# This framework has not been reproduced for the moment. 
#--------------------------------------------------------------------------------------

if Thiery_2021:
    
    sys.path.append(os.path.abspath(scripts_dir+"/figures/thiery_2021"))

    #Configuration of the plots#

    plot_fig1 = 0 # 0: do not plot figure 1 of paper
                  # 1: plot figure 1 of paper
    plot_fig2 = 0 # 0: do not plot figure 2 of paper
                  # 1: plot figure 2 of paper

    print("-------------------------------------------------------")
    print("Start plot_ms_and_si framework from Thiery et al.(2021)")
    print("-------------------------------------------------------")

    from plot_ms_and_si import *
    
    if plot_fig1 :

        print("Performing Plot f1 of Thiery et al.(2021)")
        plot_fig1()
    

#%%------------------------------------------------------------------------------------
# Framework to plots all figures associated to Grant et al.(2025) 
#--------------------------------------------------------------------------------------

if Grant_2025:

    sys.path.append(os.path.abspath(scripts_dir+"/figures/grant_2025"))

    #Configuration of the plots#

    plot_ms = True # Plots used in the main manuscript of Grant et al.(2025)
    plot_si = True # Plots used in the supplmentary materials of Grant et al.(2025)

    if plot_ms:

        print("--------------------------------------------------")
        print("Start plot_ms framework from Grant et al.(2025)")
        print("--------------------------------------------------")

        from plot_ms import *

        # f1 of ms, conceptual figure of city grid cell
        print("Performing Plot f1 of Grant et al.(2025)")
        plot_conceptual(
            da_cohort_size,
            countries_mask,
            countries_regions,
            d_isimip_meta,
            flags,
            df_life_expectancy_5,
        )

        # f2 of ms, combined heatwave plot
        print("Performing Plot f2 of Grant et al.(2025)")
        plot_combined_piechart(
            df_GMT_strj,
            ds_pf_gs,
            da_gs_popdenom,
            gdf_country_borders,
            sims_per_step,
            flags,
            df_countries,
        )
        
        # # f2 alternative with both absolute pops below box plots and pie charts
        # plot_combined_population_piechart(
        #     df_GMT_strj,
        #     ds_pf_gs,
        #     da_gs_popdenom,
        #     gdf_country_borders,
        #     sims_per_step,
        #     flags,
        #     df_countries,
        # )    
        
        # f2 alternative with absolute pops below box plots and no pie charts
        # further, returning robinson boundaries for use in pyramid plot maps for consistent map extents (that exclude antarctica)
        print("Performing Plot f2_alt of Grant et al.(2025)")
        gdf_robinson_bounds = plot_combined_population(
            df_GMT_strj,
            ds_pf_gs,
            da_gs_popdenom,
            gdf_country_borders,
            sims_per_step,
            flags,
            df_countries,
        )        

        # f3 of heatmaps across all hazards
        print("Performing Plot f3 of Grant et al.(2025)")
        plot_heatmaps_allhazards(
            df_GMT_strj,
            da_gs_popdenom,
            flags,
        )

        # # f4 of emergence union plot for hazards between 1960 and 2020 in a 2.7 degree world
        # plot_emergence_union(
        #     grid_area,
        #     da_emergence_mean,
        # )

        # # f4 alternative for hexagons and multiple thresholds
        # plot_hexagon_multithreshold(
        #     d_global_emergence,
        # )    

        # # f4 pyramid plotting
        # # first set up quantiles for plotting
        # pyramid_setup(
        #     flags,
        #     ds_gdp,
        #     ds_grdi,
        #     da_cohort_size_1960_2020,
        #     ds_vulnerability,
        # )
        # # then run plots
        # for vln_type in ('gdp','grdi'):
        #     pyramid_plot(
        #         flags,
        #         df_GMT_strj,
        #         vln_type,
        #     )

    if plot_si:

        print("--------------------------------------------------")
        print("Start plot_si framework from Grant et al.(2025)")
        print("--------------------------------------------------")

        from plot_si import *

        # heatmaps but with simulations limited to common sims (to avoid dry GCM jumps)
        print("Performing Plot sf1 of Grant et al.(2025)")
        plot_sf1_heatmaps_allhazards(
            df_GMT_strj,
            da_gs_popdenom,
            flags,
        )    
        
        # pf box plots for 1.5, 2.5 and 3.5 degree world across birth years
        print("Performing Plot sf2 of Grant et al.(2025)")
        plot_sf2_boxplots_allhazards(
            da_gs_popdenom,
            df_GMT_strj,
            flags,
        )      
        
        # pf time series for 2.7 degree world across birth years
        print("Performing Plot sf3 of Grant et al.(2025)")
        plot_sf3_pf_by_tseries_allhazards(
            flags,
            df_GMT_strj,
            da_gs_popdenom,
        )          
        
        # pf maps for 1..5, 2.5, 3.5 for all hazards
        print("Performing Plot sf4 of Grant et al.(2025)")
        plot_sf4_pf_maps_allhazards(
            da_gs_popdenom,
            gdf_country_borders,
            flags,
        )        
        
        # emergence fraction plot for hazards between 1960 and 2020 in a 2.7 degree world
        # ! Needs to flags['global_avg_emergence'] to be activated to produce this plot 
        if flags['global_avg_emergence']:
            print("Performing Plot sf5 of Grant et al.(2025)")
            plot_sf5_emergence_fracs(
                grid_area,
                ds_emergence_mean,
            )        
        else: 
            print("Can not produce Plot sf5 of Grant et al.(2025) because flags['global_avg_emergence']=0")
        
        # plot locations where exposure occurs at all in our dataset
        print("Performing Plot sf6 of Grant et al.(2025)")
        print("Plot sf6 can not be produced. Data are missing.")
        # plot_sf6_exposure_locations(
        #     grid_area,
        #     countries_mask,
        #     flags,
        # )        
        
        # plot heatmaps of pf for country level emergence
        print("Performing Plot sf7 of Grant et al.(2025)")
        plot_sf7_heatmaps_allhazards_countryemergence(
            df_GMT_strj,
            flags,
        )     
        
        # plot gmt time series for projections (rcp) and for which we map projections onto (ar6)
        print("Performing Plot sf8 of Grant et al.(2025)")
        plot_sf8_gmt_pathways(
            df_GMT_strj,
            d_isimip_meta,
        )    
            

        # pf time series for 2020 birth year across GMTs
        print("Performing Plot sfa of Grant et al.(2025)")
        print("Plot sfa can not be produced. Errors in the original script.")
        # plot_pf_gmt_tseries_allhazards(
        #     df_GMT_strj,
        #     da_gs_popdenom,
        #     flags,
        # )
        
        # plot tseries box plots for 1.5, 2.5 and 3.5 when denominator contrained by exposure extent
        print("Performing Plot sfb of Grant et al.(2025)")
        print("Plot sfb can not be produced. Data are missing.")
        # plot_geoconstrained_boxplots(
        #     flags,
        # )    
        
        # plot pie charts of all hazards
        print("Performing Plot sfc of Grant et al.(2025)")
        if flags['extr']=='all':
            plot_allhazards_piecharts(
            da_gs_popdenom,
            df_countries,
            flags,
        )
        else:
            print("Plot sfc can not be produced. Data are missing.")
        
        
        # plot cohort sizes in stacked bar chart
        print("Performing Plot sfd of Grant et al.(2025)")
        print("Plot sfb can not be produced. Data are missing.")
        # plot_cohort_sizes(
        #     df_countries,
        #     da_gs_popdenom,
        # )    
        
        # plot hexagon landfracs (will change to only show landfracs for SI)
        print("Performing Plot sfe of Grant et al.(2025)")
        print("Plot sfe can not be produced. Data are missing.")
        # plot_hexagon_landfrac(
        #     d_global_emergence,
        # )    
        
        # plot heatmaps of delta CF between main text f3 (heatwavedarea panel) and 
        print("Performing Plot sff of Grant et al.(2025)")
        print("Plot sfa can not be produced. Errors in the original script.")
        # plot_life_expectancy_testing(
        #     df_GMT_strj,
        #     GMT_indices_plot,
        #     da_gs_popdenom,
        #     flags,
        # )    

#%%------------------------------------------------------------------------------------
# Framework to plots all figures associated to Laridon et al.(2025) 
#--------------------------------------------------------------------------------------

if Laridon_2025:

    sys.path.append(os.path.abspath(scripts_dir+"/figures/laridon_2025"))

    #Configuration of the plots#

    plot_fig1 = 0 # 0: do not plot figure 1 of paper
                  # 1: plot figure 1 of paper

    print("--------------------------------------------------------")
    print("Start plot_ms_and_si framework from Laridon et al.(2025)")
    print("--------------------------------------------------------")

    from plot_ms_and_si import *
    
    if plot_fig1 :

        print("Performing Plot f1 of Laridon et al.(2021)")
        plot_fig1()

#%%------------------------------------------------------------------------------------
# Framework to plots not configured 
#--------------------------------------------------------------------------------------

if not(Thiery_2021 or Grant_2025 or Laridon_2025):
    print("No pre-defined Plots Framework")

