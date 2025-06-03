# ------------------------------------------#
# Subscripts to execute the plots functions #
# ------------------------------------------#

#%%-------------------------------------------------------------- #
# Libraries                                                       #
# --------------------------------------------------------------- #

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

#%%-----------------------------------------------------------------------#
# Framework to plots all figures associated to Thiery et al.(2021)        #
# This framework has not been reproduced for the moment.                  #
#-------------------------------------------------------------------------#
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
    

#%%-----------------------------------------------------------------------#
# Framework to plots all figures associated to Grant et al.(2025)         #
#-------------------------------------------------------------------------#

if Grant_2025:

    sys.path.append(os.path.abspath(scripts_dir+"/figures/grant_2025"))

    #Configuration of the plots#

    plot_ms = 1                     # Plots used in the main manuscript of Grant et al.(2025)

    plot_si = 0                     # Plots used in the supplmentary materials of Grant et al.(2025)

    if plot_ms:

        grant2025_fig1 = False
        grant2025_fig2 = False
        grant2025_fig2_alt = False
        grant2025_fig2_alt_mod = True
        grant2025_fig3 = False

        print("--------------------------------------------------")
        print("Start plot_ms framework from Grant et al.(2025)")
        print("--------------------------------------------------")

        from plot_ms import *

        if grant2025_fig1:

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

        if grant2025_fig2:

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
        
        if grant2025_fig2_alt:

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
        
        if grant2025_fig2_alt_mod:

            # f2 alternative with absolute pops below box plots and no pie charts
            # further, returning robinson boundaries for use in pyramid plot maps for consistent map extents (that exclude antarctica)
            print("Performing Plot f2_alt update of Grant et al.(2025)")
            gdf_robinson_bounds = plot_combined_population_update(
                df_GMT_strj,
                ds_pf_gs,
                da_gs_popdenom,
                gdf_country_borders,
                sims_per_step,
                flags,
                df_countries,
            )  

        if grant2025_fig3:

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
    
#%%-----------------------------------------------------------------------#
# Framework to plots all figures associated to Assessment                 #
#-------------------------------------------------------------------------#

# Assessment based on Grant et al.(2025) version 
if Grant_2025: 

    Norway_BiCC2 = False             # Plots used for the BiCC2 Norway law suit assessment
    Website = False                   # Plots used for the website assessment
    
    if Norway_BiCC2:

        sys.path.append(os.path.abspath(scripts_dir+"/figures/assessment"))

        from plot_assessment import *

        Norway_BiCC2_fig1 = 0      # 0: do not plot figure 1 for BiCC2 Norway law suit
                                   # 1: plot figure 1 for BiCC2 Norway law suit
    

        #--- WORK IN PROGRESS not working ---#
        Norway_BiCC2_fig2 = 0      # 0: do not plot figure 2 for BiCC2 Norway law suit
                                   # 1: plot figure 2 for BiCC2 Norway law suit
        #------------------------------------#

        if Norway_BiCC2_fig1:

            print("Performing Plot f1 BiCC2 Norway Assessmennt")

            plot_Norway_BiCC2_fig1(
                da_cohort_size,
                countries_mask,
                countries_regions,
                d_isimip_meta,
                flags,
                df_life_expectancy_5,
            )
        
        if Norway_BiCC2_fig2:

            print("Performing Plot f2 BiCC2 Norway Assessmennt")

            plot_Norway_BiCC2_fig2(
                df_GMT_strj,
                ds_pf_gs,
                da_gs_popdenom,
                gdf_country_borders,
                sims_per_step,
                flags,
                df_countries,
            )


    if Website:

        sys.path.append(os.path.abspath(scripts_dir+"/figures/assessment"))

        from plot_assessment import *

        website_fig1 = 1        # 0: do not plot figure 1 for BiCC2 Norway law suit
                                # 1: plot figure 1 for BiCC2 Norway law suit

        if website_fig1:

            print("Performing Plot f1 website assessmennt")

            plot_website_fig1(
                da_cohort_size,
                countries_mask,
                countries_regions,
                d_isimip_meta,
                flags,
                df_life_expectancy_5,
            )


#%%-----------------------------------------------------------------------#
# Framework to plots all figures associated to Laridon et al.(2025)       #
#-------------------------------------------------------------------------#

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

#%%-----------------------------------------------------------------------#
# Framework to plots all figures associated to Source2Suffering           #
#-------------------------------------------------------------------------#
if Source2Suffering:
    
    sys.path.append(os.path.abspath(scripts_dir+"/figures/source2suffering"))

    #-----------------------------------------------------------------------------------------#
    #                        Configuration of the plots to execute                            #
    #-----------------------------------------------------------------------------------------#

    # --------------- Plots for Development of the S2S Framework ----------------- #

    # Fig1 - LE per country per run for strj # 
    
    plot_fig1 = 0       # 0: do not plot figure 1
                        # 1: plot figure 1 
    
    # Fig 2 - LE trends per AR6 regions per run for strj #

    plot_fig2 = 0       # 0: do not plot figure 2 
                        # 1: plot figure 2 
    
    # Fig 3 - LE per region per run for strj #

    plot_fig3 = 0       # 0: do not plot figure 3 
                        # 1: plot figure 3

    # Fig 4 - LE MMM per country per run for strj #

    plot_fig4 = 0       # 0: do not plot figure 4 
                        # 1: plot figure 4 

    # Fig 5 - LE MMM per region per run for strj #

    plot_fig5 = 0       # 0: do not plot figure 5
                        # 1: plot figure 5 

    # Fig 6 - LE MMM and EMF for all countries and regions for 1.5/2.5/3.5°C #

    plot_fig6 = 1       # 0: do not plot figure 6
                        # 1: plot figure 6 

    # Fig 7 - BE LE MMM for all regions  #

    plot_fig7 = 1       # 0: do not plot figure 7  
                        # 1: plot figure 7 

    # Fig 8 - BE LE MMM for all countries  #

    plot_fig8 = 0       # 0: do not plot figure 8 
                        # 1: plot figure 8 

    # Fig 14 - LFE MMM for the world region  #

    plot_fig14 = 0      # 0: do not plot figure 14 
                        # 1: plot figure 14

    # Fig 15 - LFE MMM for all countries  #

    plot_fig15 = 0      # 0: do not plot figure 15 
                        # 1: plot figure 15

    # --------------- Plots for Assessment Source2Suffering - Neptun Deep -------- #

    # Fig 9 - Number of new exposed children due to specific GHG emission and comparison with flags['gmt] and Thiery et al.(2021)  #

    plot_fig9 = 0       # 0: do not plot figure 9 
                        # 1: plot figure 9 
    plot_fig10 = 0      # 0: do not plot figure 10  
                        # 1: plot figure 10  (only works if plot_fig9=1)
    plot_fig16 = 0      # 0: do not plot figure 16 
                        # 1: plot figure 16

    # -------------------- Plots for Assessment - SPARCCLE-STS ------------------- #

    # Fig 11 - Comparison of the different GMT pathways between GMT_15, GMT_20, GMT_NDC, GMT_OS and GMT_noOS  #

    plot_fig11 = 0      # 0: do not plot figure 11 
                        # 1: plot figure 11 

    # Fig 12 - Comparison of the different GMT pathways between GMT_OS, GMT_noOS and GMT_STS  #

    plot_fig12 = 0      # 0: do not plot figure 12 
                        # 1: plot figure 12 

    # Fig 13 - LE MMM for all regions under STS_ModAct and STS_Ren  #

    plot_fig13 = 0      # 0: do not plot figure 13 
                        # 1: plot figure 13
                        
    # Fig 17 - LFE MMM for all regions under STS_ModAct and STS_Ren  #

    plot_fig17 = 0      # 0: do not plot figure 17 
                        # 1: plot figure 17

    
    #-----------------------------------------------------------------------------------------#
    #                                        Init                                             #
    #-----------------------------------------------------------------------------------------#    
    print(" -------------------------------------------------------")
    print("|        Start plots for Source2Suffering Project       |")
    print(" -------------------------------------------------------")

    adr_plot_source2suffering = scripts_dir+'/figures/source2suffering/plot_source2suffering.py'
    with open(adr_plot_source2suffering) as f:
        exec(f.read(), globals())

    # Definition of the hazards 
        hazards = [
            "burntarea",
            "cropfailedarea",
            "driedarea",
            "floodedarea",
            "heatwavedarea",
            "tropicalcyclonedarea"
        ]

    # -------------------------------------------------------------------------------------------------------- #
    #                                          Plots for Development                                           # 
    # -------------------------------------------------------------------------------------------------------- #
    
    if plot_fig1:

        print("\n----- Performing Plot f1 for Development - LE per country per run for strj -----")

        for extr in hazards:

            print("\nPerforming Plot f1 for {}".format(extr))

            with open(data_dir+'{}/{}/ds_le_percountry_perrun_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                dataset_le = pk.load(f)

            #if flags['extr']=='heatwavedarea':

            # North America #

            plot_dev_fig1(flags, extr, ds=dataset_le, country='United States', run=6, GMT=0)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='United States', run=6, GMT=10)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='United States', run=6, GMT=20)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='Canada', run=6, GMT=0)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='Canada', run=6, GMT=10)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='Canada', run=6, GMT=20)

            # South Asia #

            plot_dev_fig1(flags, extr, ds=dataset_le, country='India', run=6, GMT=0)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='India', run=6, GMT=10)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='India', run=6, GMT=20)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='Afghanistan', run=6, GMT=0)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='Afghanistan', run=6, GMT=10)
            plot_dev_fig1(flags, extr, ds=dataset_le, country='Afghanistan', run=6, GMT=20)

    
    if plot_fig2:
    
        print("\n----- Performing Plot f2 for Development - LE trends per AR6 regions per run for strj -----")

        if flags['extr']=='heatwavedarea':

            print("\nPerforming Plot f2 for {}".format(extr))

            with open(data_dir+'{}/{}/lifetime_exposure_trends_regions.pkl'.format(flags['version'],flags['extr']), 'rb') as f:
                ds_le_trends_regions = pk.load(f)

            plot_dev_fig2(flags, ds=ds_le_trends_regions, var_name='exposure_trend_ar6', coords_dict={'run':6,'GMT':20,'region':17})
    
    if plot_fig3:
        
        print("\n----- Performing Plot f3 for Development - LE per region per run for strj -----")

        for extr in hazards:

            print("\nPerforming Plot f3 for {}".format(extr))

            with open(data_dir+'{}/{}/ds_le_perregion_perrun_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_regions = pk.load(f)

            #if flags['extr']=='heatwavedarea':

            # North America #

            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=7, run=6, GMT=0)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=7, run=6, GMT=10)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=7, run=6, GMT=20)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=7, run=11, GMT=0)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=7, run=11, GMT=10)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=7, run=11, GMT=20)
            
            # South Asia #

            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=8, run=6, GMT=0)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=8, run=6, GMT=10)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=8, run=6, GMT=20)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=8, run=11, GMT=0)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=8, run=11, GMT=10)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=8, run=11, GMT=20)

            # World #

            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=11, run=6, GMT=0)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=11, run=6, GMT=10)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=11, run=6, GMT=20)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=11, run=11, GMT=0)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=11, run=11, GMT=10)
            plot_dev_fig3(ds_regions, flags, extr, ds=ds_le_regions, region=11, run=11, GMT=20)
    
    if plot_fig4:
        
        print("\n----- Performing Plot f4 for Development - LE MMM per country per run for strj -----")

        for extr in hazards:

            print("\nPerforming Plot f4 for {}".format(extr))

            with open(data_dir+'{}/{}/ds_le_percountry_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_mmm_country = pk.load(f)

            # North America #

            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='United States', GMT=0)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='United States', GMT=10)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='United States', GMT=20)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='Canada', GMT=0)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='Canada', GMT=10)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='Canada', GMT=20)
            
            # South Asia #

            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='India', GMT=0)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='India', GMT=10)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='India', GMT=20)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='Afghanistan', GMT=0)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='Afghanistan', GMT=10)
            plot_dev_fig4(flags, extr, ds=ds_le_mmm_country, country='Afghanistan', GMT=20)
    
    if plot_fig5:
        
        print("\n----- Performing Plot f5 for Development - LE MMM per region per run for strj -----")

        for extr in hazards:

            print("\nPerforming Plot f5 for {}".format(extr))

            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_mmm_regions = pk.load(f)

            # North America #

            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=7, GMT=0)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=7, GMT=10)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=7, GMT=20)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=7, GMT=0)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=7, GMT=10)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=7, GMT=20)
            
            # South Asia #

            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=8, GMT=0)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=8, GMT=10)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=8, GMT=20)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=8, GMT=0)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=8, GMT=10)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=8, GMT=20)

            # World #

            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=11, GMT=0)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=11, GMT=10)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=11, GMT=20)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=11, GMT=0)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=11, GMT=10)
            plot_dev_fig5(ds_regions, flags, extr, ds=ds_le_mmm_regions, region=11, GMT=20)

    if plot_fig6:

        print("\n----- Performing Plot f6 for Story Lines for myclimatefuture - LE MMM & EMF for all countries and regions for 1.5/2.5/3.5°C -----")

        for extr in hazards:

            #------------------------------------ Regions ----------------------------------------#

            print("\nPerforming Plot f6 for {} for regions\n".format(extr))

            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format('pickles_sandbox',extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_perregion = pk.load(f)
            
            # with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
            #     ds_le_perregion = pk.load(f)

            for region_ind in ds_le_perregion['mmm_BE']['region'].values:
        
                plot_dev_fig6_region(ds_regions, flags, extr, ds=ds_le_perregion, region=region_ind, EMF=False)
            
            # with open(data_dir+'{}/{}/ds_EMF_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt']), 'rb') as f:
            #     ds_EMF_mmm_region = pk.load(f)

            # for region_ind in ds_le_perregion['mmm_BE']['region'].values:
        
            #     plot_dev_fig6_region(ds_regions, flags, extr, ds=ds_EMF_mmm_region, region=region_ind, EMF=True)

            #------------------------------------ Country ----------------------------------------#

            # print("\nPerforming Plot f6 for {} for countries\n".format(extr))

            # with open(data_dir+'{}/{}/ds_le_percountry_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
            #     ds_le_percountry = pk.load(f)

            # for country_name in ds_le_percountry['country'].values:
        
            #     plot_dev_fig6_country(flags, extr, ds=ds_le_percountry, country=country_name, EMF=False)

            # with open(data_dir+'{}/{}/ds_EMF_percountry_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt']), 'rb') as f:
            #     ds_EMF_mmm_country = pk.load(f)

            # for country_name in ds_le_percountry['country'].values:
        
            #     plot_dev_fig6_country(flags, extr, ds=ds_EMF_mmm_country, country=country_name, EMF=True)

    
    if plot_fig7:

        print("\n-----Performing Plot f7 for Development - BE LE MMM for all regions -----")

        for extr in hazards:

            print("\nPerforming Plot f7 for {}\n".format(extr))

            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format('pickles_sandbox',extr,flags['gmt'],flags['rm']), 'rb') as f:
                    ds_le_perregion = pk.load(f)

            # with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
            #         ds_le_perregion = pk.load(f)

            for i in ds_regions['region'].values:

                plot_dev_fig7(ds_regions, df_GMT_strj, ds_le_perregion, extr, flags, region_ind=i)
    
    if plot_fig8:
    
        print("\n----- Performing Plot f8 for Development - BE LE MMM for all countries -----")

        for extr in hazards:

            print("\nPerforming Plot f8 for {}\n".format(extr))

            with open(data_dir+'{}/{}/ds_le_percountry_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_percountry = pk.load(f)

            for country_name in ds_le_percountry['country'].values:
    
                plot_dev_fig8(df_GMT_strj, ds_le_percountry, flags, extr, country=country_name)

    if plot_fig14:
        
        print("\n----- Performing Plot f14 for Development : LFE MMM for the World region -----")

        if flags['gmt'] =='ar6_new':
    
            print("\nPlot not performed because routine with the stylized trajectories not yet settle for the LFE with flags['gmt'] = ar6_new")

        if flags['gmt']=='original':

            for extr in hazards:

                print("\nPerforming Plot f14 for {}".format(extr))

                with open(data_dir+'{}/{}/ds_lfe_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                    ds_lfe_perregion = pk.load(f)

                plot_dev_fig14_regions(flags, extr=extr, ds_lfe_perregion=ds_lfe_perregion, ind_region=11)
        

    if plot_fig15:
        
        print("\n----- Performing Plot f15 for Development : LFE MMM for all countries -----")

        if flags['gmt'] =='ar6_new':

            print("\nPlot not performed for because routine with the stylized trajectories not yet settle for the LFE with flags['gmt'] = ar6_new")

        if flags['gmt']=='original':

            for extr in hazards:

                print("\nPerforming Plot f15 for {}\n".format(extr))

                with open(data_dir+'{}/{}/ds_lfe_percountry_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                    ds_lfe_percountry = pk.load(f)
            
                for country_name in ds_le_percountry['country'].values:

                    plot_dev_fig15_country(flags, extr=extr, ds_lfe_percountry=ds_lfe_percountry, country_name=country_name)

    # -------------------------------------------------------------------------------------------------------- #
    #                           Plots for Assessment Source2Suffering- Neptun Deep                             # 
    # -------------------------------------------------------------------------------------------------------- #

    if plot_fig9:
        
        print("\n----- Performing Plot f9 for Assessment - Neptun Deep : Number of new exposed children to hazard due to specific GHG emission -----")

        # Load the values computed by the Source2Suffering framework #

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],flags['gmt'],flags['rm']), 'rb') as f:
            ds_S2S_NeptunDeep = pk.load(f)

        # Load the values computed by the Wim Thiery Expert Report #

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_WT_NeptunDeep.pkl'.format(flags['version']), 'rb') as f:
            ds_WT_NeptunDeep = pk.load(f)

        # Plots 

        for extr in ds_S2S_NeptunDeep.coords["hazard"]:

            ds_S2S_NeptunDeep_sel = ds_S2S_NeptunDeep['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)                          
            ds_WT_NeptunDeep_sel = ds_WT_NeptunDeep['wt_valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr)

            plot_dev_fig9(ds_S2S_NeptunDeep_sel,ds_WT_NeptunDeep_sel,birth_cohort_int,extr,flags)

    if plot_fig10:

        print("\n----- Performing Plot f10 for Assessment - Neptun Deep : Number of new exposed children to hazard comparison with flags['gmt] and Thiery et al.(2021) -----")

        # Load the values computed by the Source2Suffering framework #

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],'ar6_new',flags['rm']), 'rb') as f:
            ds_S2S_NeptunDeep_gmt_ar6_new_flags_rm = pk.load(f)

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],'original',flags['rm']), 'rb') as f:
            ds_S2S_NeptunDeep_gmt_original_flags_rm = pk.load(f)

        # Load the values computed by the Wim Thiery Expert Report #

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_WT_NeptunDeep.pkl'.format(flags['version']), 'rb') as f:
            ds_WT_NeptunDeep = pk.load(f)
        
        for extr in ds_WT_NeptunDeep.coords["hazard"]:
    
            ds_S2S_NeptunDeep_gmt_ar6_new_flags_rm_sel = ds_S2S_NeptunDeep_gmt_ar6_new_flags_rm['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)
            ds_S2S_NeptunDeep_gmt_original_flags_rm_sel = ds_S2S_NeptunDeep_gmt_original_flags_rm['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)                       
            ds_WT_NeptunDeep_sel = ds_WT_NeptunDeep['wt_valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr)

            plot_dev_fig10(ds_S2S_NeptunDeep_gmt_ar6_new_flags_rm_sel, ds_S2S_NeptunDeep_gmt_original_flags_rm_sel, ds_WT_NeptunDeep_sel, birth_cohort_int, extr, flags)
    
    if plot_fig16:

        print("\n----- Performing Plot f16 for Assessment - Neptun Deep : Number of new exposed children to hazard : Comparison between the 'rm' and 'gmt' configurations -----")

        # Load the values computed by the Wim Thiery Expert Report #

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_WT_NeptunDeep.pkl'.format(flags['version']), 'rb') as f:
            ds_WT_NeptunDeep = pk.load(f)
        
        # Load the values computed by the Source2Suffering framework #

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],'ar6_new','rm'), 'rb') as f:
            ds_S2S_NeptunDeep_gmt_ar6_new_rm = pk.load(f)
        
        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],'ar6_new','no_rm'), 'rb') as f:
            ds_S2S_NeptunDeep_gmt_ar6_new_no_rm = pk.load(f)

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],'original','rm'), 'rb') as f:
            ds_S2S_NeptunDeep_gmt_original_rm = pk.load(f)

        with open(data_dir+'{}/source2suffering/assessment/Neptun_Deep/ds_S2S_NeptunDeep_gmt_{}_{}.pkl'.format(flags['version'],'original','no_rm'), 'rb') as f:
            ds_S2S_NeptunDeep_gmt_original_no_rm = pk.load(f)

        year_start_as = 2010
        year_end_as = 2020

        birth_cohort_int = list(range(year_start_as, year_end_as + 1))
        
        for extr in ds_S2S_NeptunDeep.coords["hazard"]:

            ds_WT_NeptunDeep_sel = ds_WT_NeptunDeep['wt_valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr)

            ds_S2S_NeptunDeep_gmt_ar6_new_rm_sel = ds_S2S_NeptunDeep_gmt_ar6_new_rm['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)
            ds_S2S_NeptunDeep_gmt_ar6_new_no_rm_sel = ds_S2S_NeptunDeep_gmt_ar6_new_no_rm['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)
            ds_S2S_NeptunDeep_gmt_original_rm_sel = ds_S2S_NeptunDeep_gmt_original_rm['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)
            ds_S2S_NeptunDeep_gmt_original_no_rm_sel = ds_S2S_NeptunDeep_gmt_original_no_rm['valc_nr_children_facing_extra_hazard_NeptunDeep'].sel(hazard=extr,birth_year=birth_cohort_int)

            plot_dev_fig16(ds_WT_NeptunDeep_sel,
             ds_S2S_NeptunDeep_gmt_ar6_new_rm_sel, 
             ds_S2S_NeptunDeep_gmt_ar6_new_no_rm_sel,
             ds_S2S_NeptunDeep_gmt_original_rm_sel,
             ds_S2S_NeptunDeep_gmt_original_no_rm_sel,
             birth_cohort_int,extr)

    # -------------------------------------------------------------------------------------------------------- #
    #                                    Plots for Assessment - SPARCCLE STS                                   # 
    # -------------------------------------------------------------------------------------------------------- #
    
    if plot_fig11:

        print("\n----- Performing Plot f11 for Assessment - SPARCCLE-STS : Comparison of the different GMT pathways between GMT_15, GMT_20, GMT_NDC, GMT_OS and GMT_noOS -----")

        plot_dev_fig11(df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_OS, df_GMT_noOS)

    if plot_fig12:

        print("\n----- Performing Plot f12 for Assessment - SPARCCLE-STS : Comparison of the different GMT pathways between GMT_OS, GMT_noOS and GMT_STS -----")

        plot_dev_fig12(df_GMT_OS, df_GMT_noOS, ds_GMT_STS)

    if plot_fig13:
    
        print("\n----- Performing Plot f13 for Assessment - SPARCCLE-STS : LE MMM for all regions under STS_ModAct and STS_Ren -----")

        # Definition of the coordinates 
        hazards = [
            "burntarea",
            "cropfailedarea",
            "driedarea",
            "floodedarea",
            "heatwavedarea",
            "tropicalcyclonedarea"
        ]

        for extr in hazards:

            # Loading the pickles 
                    
            with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_le_perregion = pk.load(f)

            for region_ind in ds_le_perregion['mmm_STS_ModAct']['region'].values:

                # Produce the figures of the absolute value of the lifetime exposure per birth cohort and save the results in a .nc file #
                
                plot_dev_fig13_regions(ds_regions, extr, flags, ds_le=ds_le_perregion, region_ind=region_ind, EMF=False)
            
            # Save the results of interest for SPARCCLE STS in a .nc file

            ds_subset_ModAct = ds_le_perregion[["mmm_STS_ModAct", "std_STS_ModAct", "median_STS_ModAct", "lqntl_STS_ModAct", "uqntl_STS_ModAct"]]
            ds_subset_Ren = ds_le_perregion[["mmm_STS_Ren", "std_STS_Ren", "median_STS_Ren", "lqntl_STS_Ren", "uqntl_STS_Ren"]]

            ds_subset_ModAct.attrs['README'] = """
            This dataset contains the values of the lifetime exposure per birth year (from 1960 to 2020) for the STS ModAct scenario for 12 different regions.

            The index of the 'region' coordinate is the following:
            0: East Asia & Pacific
            1: Europe & Central Asia
            2: High income
            3: Latin America & Caribbean
            4: Low income
            5: Lower middle income
            6: Middle East & North Africa
            7: North America
            8: South Asia
            9: Sub-Saharan Africa
            10: Upper middle income
            11: World

            For the variables 'mmm' = multi model mean, 'std' = standard deviation, 'median' = median (p0.5), 'lqntl' = lower quartile (p0.25), 'uqntl' = upper quartile (p0.75)

            Contact author : Amaury.Laridon@vub.be
            """

            ds_subset_Ren.attrs['README'] = """
            This dataset contains the values of the lifetime exposure per birth year (from 1960 to 2020) for the STS Ren scenario for 12 different regions.

            The index of the 'region' coordinate is the following:
            0: East Asia & Pacific
            1: Europe & Central Asia
            2: High income
            3: Latin America & Caribbean
            4: Low income
            5: Lower middle income
            6: Middle East & North Africa
            7: North America
            8: South Asia
            9: Sub-Saharan Africa
            10: Upper middle income
            11: World

            For the variables 'mmm' = multi model mean, 'std' = standard deviation, 'median' = median (p0.5), 'lqntl' = lower quartile (p0.25), 'uqntl' = upper quartile (p0.75)

            Contact author : Amaury.Laridon@vub.be
            """

            ds_subset_ModAct.to_netcdf(scripts_dir + '/output/assessment/SPARCCLE_STS/STS_ModAct_{}_lifetime_exposure_perregion_{}.nc'.format(extr,flags['rm']))
            ds_subset_Ren.to_netcdf(scripts_dir + '/output/assessment/SPARCCLE_STS/STS_Ren_{}_lifetime_exposure_perregion_{}.nc'.format(extr,flags['rm']))

    if plot_fig17:
        
        print("\n----- Performing Plot f17 for Assessment - SPARCCLE-STS : LFE MMM for all regions under STS_ModAct and STS_Ren -----")

        # Definition of the coordinates 
        hazards = [
            "burntarea",
            "cropfailedarea",
            "driedarea",
            "floodedarea",
            "heatwavedarea",
            "tropicalcyclonedarea"
        ]

        for extr in hazards:

            # Loading the pickles 
                    
            with open(data_dir+'{}/{}/ds_lfe_perregion_gmt_{}_{}.pkl'.format(flags['version'],extr,flags['gmt'],flags['rm']), 'rb') as f:
                ds_lfe_perregion = pk.load(f)

            for region_ind in ds_lfe_perregion['mmm_STS_ModAct']['region'].values:

                # Produce the figures of the absolute value of the lifetime exposure per birth cohort and save the results in a .nc file #
                
                plot_dev_fig17_regions(ds_regions, extr, flags, ds_lfe=ds_lfe_perregion, region_ind=region_ind, EMF=False)
            
            # Save the results of interest for SPARCCLE STS in a .nc file

            ds_subset_ModAct = ds_lfe_perregion[["mmm_STS_ModAct", "std_STS_ModAct", "median_STS_ModAct","mmm_STS_ModAct_sm", "std_STS_ModAct_sm", "median_STS_ModAct_sm", "lqntl_STS_ModAct", "uqntl_STS_ModAct"]]
            ds_subset_Ren = ds_lfe_perregion[["mmm_STS_Ren", "std_STS_Ren", "median_STS_Ren", "mmm_STS_Ren_sm", "std_STS_Ren_sm", "median_STS_Ren_sm", "lqntl_STS_Ren", "uqntl_STS_Ren"]]

            ds_subset_ModAct.attrs['README'] = """
            This dataset contains the values of the Land Fraction Exposed (LFE) per birth year (from 1960 to 2020) for the STS ModAct scenario for 12 different regions.

            The index of the 'region' coordinate is the following:
            0: East Asia & Pacific
            1: Europe & Central Asia
            2: High income
            3: Latin America & Caribbean
            4: Low income
            5: Lower middle income
            6: Middle East & North Africa
            7: North America
            8: South Asia
            9: Sub-Saharan Africa
            10: Upper middle income
            11: World

            For the variables 'mmm' = multi model mean, 'std' = standard deviation, 'median' = median (p0.5), 'lqntl' = lower quartile (p0.25), 'uqntl' = upper quartile (p0.75). 
            'sm' means there is a 10 year rolling mean applied to smooth the results.

            Contact author : Amaury.Laridon@vub.be
            """

            ds_subset_Ren.attrs['README'] = """
            This dataset contains the values of the Land Fraction Exposed (LFE) per birth year (from 1960 to 2020) for the STS Ren scenario for 12 different regions.

            The index of the 'region' coordinate is the following:
            0: East Asia & Pacific
            1: Europe & Central Asia
            2: High income
            3: Latin America & Caribbean
            4: Low income
            5: Lower middle income
            6: Middle East & North Africa
            7: North America
            8: South Asia
            9: Sub-Saharan Africa
            10: Upper middle income
            11: World

            For the variables 'mmm' = multi model mean, 'std' = standard deviation, 'median' = median (p0.5), 'lqntl' = lower quartile (p0.25), 'uqntl' = upper quartile (p0.75)
            'sm' means there is a 10 year rolling mean applied to smooth the results.

            Contact author : Amaury.Laridon@vub.be
            """

            ds_subset_ModAct.to_netcdf(scripts_dir + '/output/assessment/SPARCCLE_STS/STS_ModAct_{}_landfraction_exposed_perregion_{}.nc'.format(extr,flags['rm']))
            ds_subset_Ren.to_netcdf(scripts_dir + '/output/assessment/SPARCCLE_STS/STS_Ren_{}_landfraction_exposed_perregion_{}.nc'.format(extr,flags['rm']))

    

#%%-----------------------------------------------------------------------#
# Framework to plots not configured                                       #
#-------------------------------------------------------------------------#

if not(Thiery_2021 or Grant_2025 or Laridon_2025 or Source2Suffering):
    print("No pre-defined Plots Framework outside Thiery et al.(2021), Grant et al.(2025), Laridon et al.(2025) or Source2Suffering Project.")

