# -------------------------------------------------------------------------------#
# Functions needed for external plots needed in the assessment reports           #
# -------------------------------------------------------------------------------#

#               
#%%---------------------------------------------------------------#
# Libraries                                                       #
# ----------------------------------------------------------------#

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
import cartopy.crs as ccrs
import seaborn as sns
import cartopy as cr
import cartopy.feature as feature
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()


#%%---------------------------------------------------------------#
# Plots for the BiCC2 Report                                      #
#-----------------------------------------------------------------#

#%%-----------------------------------------------------------------#
# Conceptual plot for emergence in one location                     #
# based on Fig.1 of Grant et al.(2025-)                             #
#-------------------------------------------------------------------#
   

def plot_Norway_BiCC2_fig1(
    da_cohort_size,
    countries_mask,
    countries_regions,
    d_isimip_meta,
    flags,
    df_life_expectancy_5,
):
    # get data
    cntry=countries
    city_name='Oslo'
    # concept_bys = np.arange(1960,2021,30)
    concept_bys = np.arange(1960,2021,1)
    print("Conceptual plot for "+cntry)
    da_smple_cht = da_cohort_size.sel(country=cntry) # cohort absolute sizes in sample country
    da_smple_cht_prp = da_smple_cht / da_smple_cht.sum(dim='ages') # cohort relative sizes in sample country
    da_cntry = xr.DataArray(
        np.in1d(countries_mask,countries_regions.map_keys(cntry)).reshape(countries_mask.shape),
        dims=countries_mask.dims,
        coords=countries_mask.coords,
    )
    da_cntry = da_cntry.where(da_cntry,drop=True)
    # weights for latitude (probably won't use but will use population instead)
    lat_weights = np.cos(np.deg2rad(da_cntry.lat))
    lat_weights.name = "weights"   
    # Oslo coords  
    city_lat = 59.9139
    city_lon = 10.7522  

    ds_spatial = xr.Dataset(
        data_vars={
            'cumulative_exposure': (
                ['run','GMT','birth_year','time','lat','lon'],
                np.full(
                    (len(list(d_isimip_meta.keys())),
                    len(GMT_indices_plot),
                    len(concept_bys),
                    len(year_range),
                    len(da_cntry.lat.data),
                    len(da_cntry.lon.data)),
                    fill_value=np.nan,
                ),
            ),
        },
        coords={
            'lat': ('lat', da_cntry.lat.data),
            'lon': ('lon', da_cntry.lon.data),
            'birth_year': ('birth_year', concept_bys),
            'time': ('time', year_range),
            'run': ('run', np.arange(1,len(list(d_isimip_meta.keys()))+1)),
            'GMT': ('GMT', GMT_indices_plot)
        }
    )

    # load demography pickle
    with open(data_dir+'{}/gridscale/gridscale_dmg_{}.pkl'.format(flags['version'],cntry), 'rb') as f:
        ds_dmg = pk.load(f)                  

    # loop over simulations
    for i in list(d_isimip_meta.keys()): 

        print('Loading simulation {} of {}'.format(i,len(d_isimip_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA = pk.load(f)
            
        # mask to sample country and reduce spatial extent
        da_AFA = da_AFA.where(ds_dmg['country_extent']==1,drop=True)
        
        for step in GMT_indices_plot:
            
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                
                da_AFA_step = da_AFA.reindex(
                    {'time':da_AFA['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
                ).assign_coords({'time':year_range})                     
                                    
                # simple lifetime exposure sum
                da_le = xr.concat(
                    [(da_AFA_step.loc[{'time':np.arange(by,ds_dmg['death_year'].sel(birth_year=by).item()+1)}].cumsum(dim='time') +\
                    da_AFA_step.sel(time=ds_dmg['death_year'].sel(birth_year=by).item()) *\
                    (ds_dmg['life_expectancy'].sel(birth_year=by).item() - np.floor(ds_dmg['life_expectancy'].sel(birth_year=by)).item()))\
                    for by in concept_bys],
                    dim='birth_year',
                ).assign_coords({'birth_year':concept_bys})
                
                da_le = da_le.reindex({'time':year_range})
                
                ds_spatial['cumulative_exposure'].loc[{
                    'run':i,
                    'GMT':step,
                    'birth_year':concept_bys,
                    'time':year_range,
                    'lat':ds_dmg['country_extent'].lat.data,
                    'lon':ds_dmg['country_extent'].lon.data,
                }] = da_le.loc[{
                    'birth_year':concept_bys,
                    'time':year_range,
                    'lat':ds_dmg['country_extent'].lat.data,
                    'lon':ds_dmg['country_extent'].lon.data,
                }]

    # mean for the city           
    da_test_city = ds_spatial['cumulative_exposure'].sel({'lat':city_lat,'lon':city_lon},method='nearest').mean(dim='run')
    da_test_city = da_test_city.rolling(time=5,min_periods=5).mean()

    # standard deviation for the city
    da_test_city_std = ds_spatial['cumulative_exposure'].sel({'lat':city_lat,'lon':city_lon},method='nearest').std(dim='run')
    da_test_city_std = da_test_city_std.rolling(time=5,min_periods=5).mean()

    # fill in 1st 4 years with 1s
    # first for mean
    for by in da_test_city.birth_year.data:
        for step in GMT_indices_plot:
            da_test_city.loc[{'birth_year':by,'GMT':step,'time':np.arange(by,by+5)}] = da_test_city.loc[{'birth_year':by,'GMT':step}].min(dim='time')
    # then for std        
    for by in da_test_city_std.birth_year.data:
        for step in GMT_indices_plot:
            da_test_city_std.loc[{'birth_year':by,'GMT':step,'time':np.arange(by,by+5)}] = da_test_city_std.loc[{'birth_year':by,'GMT':step}].min(dim='time')        
                
    # load PIC pickles
    with open(data_dir+'{}/{}/{}/gridscale_le_pic_{}_{}.pkl'.format(flags['version'],flags['extr'],cntry,flags['extr'],cntry), 'rb') as f:
        ds_pic = pk.load(f)   
    with open(data_dir+'{}/{}/{}/gridscale_pic_qntls_{}_{}.pkl'.format(flags['version'],flags['extr'],cntry,flags['extr'],cntry), 'rb') as f:
        ds_pic_qntl = pk.load(f)

    # plotting city lat/lon pixel doesn't give smooth kde
    df_pic_city = ds_pic['lifetime_exposure'].sel({'lat':city_lat,'lon':city_lon},method='nearest').to_dataframe().drop(columns=['lat','lon',])         
    da_pic_city_9999 = ds_pic_qntl['99.99'].sel({'lat':city_lat,'lon':city_lon},method='nearest')  

    # concept figure
    # ------------------------------------------------------------------   
    
    # plot building
    from mpl_toolkits.axes_grid1 import inset_locator as inset
    plt.rcParams['patch.linewidth'] = 0.1
    plt.rcParams['patch.edgecolor'] = 'k'
    colors = dict(zip(GMT_indices_plot,['steelblue','darkgoldenrod','darkred']))
    x=5
    y=1
    l = 0
    gmt_legend={
        GMT_indices_plot[0]:'1.5',
        GMT_indices_plot[1]:'2.5',
        GMT_indices_plot[2]:'3.5',
    }

    f,ax = plt.subplots(
        figsize=(x,y)
        )

    # ------------------------------------------------------------------   
    # 1960 time series
    

    for step in GMT_indices_plot:
        da_test_city.loc[{'birth_year':1960,'GMT':step}].plot.line(
            ax=ax,
            color=colors[step],
            linewidth=1,
        )
        # bold line for emergence
        da = da_test_city.loc[{'birth_year':1960,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        da.plot.line(
            ax=ax,
            color=colors[step],
            linewidth=3,
            zorder=4,
        )
             
    end_year=1960+np.floor(df_life_expectancy_5.loc[1960,cntry])
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    ax.set_xticks(np.arange(1960,2031,10))
    ax.set_xticklabels([1960,None,1980,None,2000,None,2020,None])
    ax.set_yticks([0,5])
    ax.set_yticklabels([None,5])     
    ax.annotate(
        'Born in 1960',
        (1965,ax.get_ylim()[-1]+2),
        xycoords=ax.transData,
        fontsize=10,
        fontweight='bold',
        rotation='horizontal',
        color='gray',
    )    
    ax.set_title(None)
    ax.annotate(
        letters[l],
        (1960,ax.get_ylim()[-1]+2),
        xycoords=ax.transData,
        fontsize=10,
        rotation='horizontal',
        color='k',
        fontweight='bold',
    )    
    l+=1      
      
    ax.set_xlim(
        1960,
        end_year,
    )
    ax.set_ylim(
        0,
        da_pic_city_9999+1,
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)    
    ax.tick_params(colors='gray')
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_color('gray')
    ax.hlines(
        y=da_pic_city_9999, 
        xmin=1960, 
        xmax=da_test_city.loc[{'birth_year':1960}].time.max()+10, 
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )

    # 1960 pdf
    ax_pdf_l = end_year+5
    ax_pdf_b = -2
    ax_pdf_w = 20
    ax_pdf_h = ax.get_ylim()[-1]+2
    ax_pdf = ax.inset_axes(
        bounds=(ax_pdf_l, ax_pdf_b, ax_pdf_w, ax_pdf_h),
        transform=ax.transData,
    )
    sns.histplot(
        data=df_pic_city.round(),
        y='lifetime_exposure',
        color='lightgrey',
        discrete = True,
        ax=ax_pdf
    )
    ax_pdf.hlines(
        y=da_pic_city_9999, 
        xmin=0, 
        xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    for step in GMT_indices_plot:
        ax_pdf.hlines(
            y=da_test_city.loc[{'birth_year':1960,'GMT':step}].max(), 
            xmin=0, 
            xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
            colors=colors[step], 
            linewidth=1, 
            linestyle='-', 
            label=gmt_legend[step], 
            zorder=2
        )
    ax_pdf.spines['right'].set_visible(False)
    ax_pdf.spines['top'].set_visible(False)      
    ax_pdf.set_ylabel(None)
    ax_pdf.set_xlabel(None)
    ax_pdf.set_ylim(-2,ax.get_ylim()[-1])
    ax_pdf.tick_params(colors='gray')
    ax_pdf.spines['left'].set_color('gray')
    ax_pdf.spines['bottom'].set_color('gray')
    ax_pdf.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10,
    )
    l+=1
        
    # ------------------------------------------------------------------       
    # 1990 time series
    ax2_l = 1990
    ax2_b = da_pic_city_9999 *2
    ax2_w = np.floor(df_life_expectancy_5.loc[1990,cntry])
    ax2_h = np.round(da_test_city.loc[{'birth_year':1990,'GMT':GMT_indices_plot[-1]}].max())
    ax2 = ax.inset_axes(
        bounds=(ax2_l, ax2_b, ax2_w, ax2_h),
        transform=ax.transData,
    )

    for step in GMT_indices_plot:
        da_test_city.loc[{'birth_year':1990,'GMT':step}].plot.line(
            ax=ax2,
            color=colors[step],
            linewidth=1,
        )
        # bold line for emergence
        da = da_test_city.loc[{'birth_year':1990,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        da.plot.line(
            ax=ax2,
            color=colors[step],
            linewidth=3,
            zorder=4,
        )    
                  
    end_year=1990+np.floor(df_life_expectancy_5.loc[1990,cntry])
    ax2.set_ylabel(None)
    ax2.set_xlabel(None)
    ax2.set_yticks([0,5,10])
    ax2.set_yticklabels([None,5,10])  
    ax2.set_xticks(np.arange(1990,2071,10))      
    ax2.set_xticklabels([None,2000,None,2020,None,2040,None,2060,None])
    ax2.set_xlim(
        1990,
        end_year,
    )
    ax2.set_ylim(
        0,
        np.round(da_test_city.loc[{'birth_year':1990,'GMT':GMT_indices_plot[-1]}].max())+1,
    )
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)  
    ax2.spines['left'].set_position(('data',1990)) 
    ax2.tick_params(colors='gray')
    ax2.spines['left'].set_color('gray')
    ax2.spines['bottom'].set_color('gray')
    
    ax2.annotate(
        'Born in 1990',
        (1995,ax2.get_ylim()[-1]),
        xycoords=ax2.transData,
        fontsize=10,
        fontweight='bold',
        rotation='horizontal',
        color='gray',
    )
    ax2.set_title(None)
    ax2.annotate(
        letters[l],
        (1990,ax2.get_ylim()[-1]),
        xycoords=ax2.transData,
        fontsize=10,
        rotation='horizontal',
        color='k',
        fontweight='bold',
    )     
    l+=1           

    # get time of first line to cross PIC thresh
    emergences = []
    for step in GMT_indices_plot:
        da = da_test_city.loc[{'birth_year':1990,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        if np.any(da.notnull()): # testing for median
            da_t = da.time.where(da == da.min()).dropna(dim='time').item()
            emergences.append(da_t)
    first_emerge = np.min(emergences)

    ax2.hlines(
        y=da_pic_city_9999, 
        xmin=first_emerge, 
        xmax=end_year, 
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )        

    # 1990 pdf
    ax2_pdf_l = end_year+5
    ax2_pdf_b = -2
    ax2_pdf_w = 20
    ax2_pdf_h = ax2.get_ylim()[-1]+2
    ax2_pdf = ax2.inset_axes(
        bounds=(ax2_pdf_l, ax2_pdf_b, ax2_pdf_w, ax2_pdf_h),
        transform=ax2.transData,
    )
    sns.histplot(
        data=df_pic_city.round(),
        y='lifetime_exposure',
        color='lightgrey',
        discrete = True,
        ax=ax2_pdf
    )
    ax2_pdf.hlines(
        y=da_pic_city_9999, 
        xmin=0, 
        xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    for step in GMT_indices_plot:
        ax2_pdf.hlines(
            y=da_test_city.loc[{'birth_year':1990,'GMT':step}].max(), 
            xmin=0, 
            xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
            colors=colors[step], 
            linewidth=1, 
            linestyle='-', 
            label=gmt_legend[step], 
            zorder=2
        )
    ax2_pdf.spines['right'].set_visible(False)
    ax2_pdf.spines['top'].set_visible(False)      
    ax2_pdf.set_ylabel(None)
    ax2_pdf.set_xlabel(None)
    ax2_pdf.set_ylim(-2,ax2.get_ylim()[-1])
    ax2_pdf.tick_params(colors='gray')
    ax2_pdf.spines['left'].set_color('gray')
    ax2_pdf.spines['bottom'].set_color('gray')
    ax2_pdf.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10,
    )
    l+=1   
    
    ax2_pdf.annotate(
        'Unprecedented\nlifetime\nexposure\nfor {} people'.format(str(int(np.round(ds_dmg['by_population_y0'].sel({'birth_year':1990,'lat':city_lat,'lon':city_lon},method='nearest').item(),-3)))),
        (1.1,0.3),
        xycoords=ax2_pdf.transAxes,
        fontsize=14,
        rotation='horizontal',
        color='gray',
        # fontweight='bold',
    )             

    # ------------------------------------------------------------------   
    # 2020 time series
    ax3_l = 2020
    ax3_b = np.round(da_test_city.loc[{'birth_year':1990,'GMT':GMT_indices_plot[-1]}].max()) * 1.5
    ax3_w = np.floor(df_life_expectancy_5.loc[2020,cntry])
    ax3_h = np.round(da_test_city.loc[{'birth_year':2020,'GMT':GMT_indices_plot[-1]}].max())
    ax3 = ax2.inset_axes(
        bounds=(ax3_l, ax3_b, ax3_w, ax3_h),
        transform=ax2.transData,
    )
    # plot mean lines
    for step in GMT_indices_plot:
        da_test_city.loc[{'birth_year':2020,'GMT':step}].plot.line(
            ax=ax3,
            color=colors[step],
            linewidth=1,
        )
        # bold line for emergence
        da = da_test_city.loc[{'birth_year':2020,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        da.plot.line(
            ax=ax3,
            color=colors[step],
            linewidth=3,
            zorder=4,
        )    

    end_year=2020+np.floor(df_life_expectancy_5.loc[2020,cntry])      
    
    ax3.set_ylabel(None)
    ax3.set_xlabel(None)
    ax3.set_yticks([0,5,10,15,20,25])
    ax3.set_yticklabels([None,5,10,15,20,25])   
    ax3.set_xticks(np.arange(2020,2101,10))      
    ax3.set_xticklabels([2020,None,2040,None,2060,None,2080,None,2100])
    ax3.set_xlim(
        2020,
        end_year,
    )
    ax3.set_ylim(
        0,
        np.round(da_test_city.loc[{'birth_year':2020,'GMT':GMT_indices_plot[-1]}].max())+1,
    )
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)  
    ax3.spines['left'].set_position(('data',2020))  
    ax3.tick_params(colors='gray')
    ax3.spines['left'].set_color('gray')
    ax3.spines['bottom'].set_color('gray')

    # get time of first line to cross PIC thresh
    emergences = []
    for step in GMT_indices_plot:
        da = da_test_city.loc[{'birth_year':2020,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        if np.any(da.notnull()):
            da_t = da.time.where(da == da.min()).dropna(dim='time').item()
            emergences.append(da_t)
    first_emerge = np.min(emergences)

    ax3.hlines(
        y=da_pic_city_9999, 
        xmin=first_emerge, 
        xmax=end_year, 
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    ax3.annotate(
        'Born in 2020',
        (2025,ax3.get_ylim()[-1]),
        xycoords=ax3.transData,
        fontsize=10,
        fontweight='bold',
        rotation='horizontal',
        color='gray',
    )
    ax3.set_title(None)
    ax3.annotate(
        letters[l],
        (2020,ax3.get_ylim()[-1]),
        xycoords=ax3.transData,
        fontsize=10,
        rotation='horizontal',
        color='k',
        fontweight='bold',
    ) 
    l+=1      

    # 2020 pdf
    ax3_pdf_l = end_year+5
    ax3_pdf_b = -2
    ax3_pdf_w = 20
    ax3_pdf_h = ax3.get_ylim()[-1]+2
    ax3_pdf = ax3.inset_axes(
        bounds=(ax3_pdf_l, ax3_pdf_b, ax3_pdf_w, ax3_pdf_h),
        transform=ax3.transData,
    )
    sns.histplot(
        data=df_pic_city.round(),
        y='lifetime_exposure',
        color='lightgrey',
        discrete = True,
        ax=ax3_pdf
    )
    ax3_pdf.hlines(
        y=da_pic_city_9999, 
        xmin=0, 
        xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    for step in GMT_indices_plot:
        ax3_pdf.hlines(
            y=da_test_city.loc[{'birth_year':2020,'GMT':step}].max(), 
            xmin=0, 
            xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
            colors=colors[step], 
            linewidth=1, 
            linestyle='-', 
            label=gmt_legend[step], 
            zorder=2
        )
    ax3_pdf.spines['right'].set_visible(False)
    ax3_pdf.spines['top'].set_visible(False)      
    ax3_pdf.set_ylabel(None)
    ax3_pdf.set_xlabel(None)
    ax3_pdf.set_ylim(-2,ax3.get_ylim()[-1])
    ax3_pdf.tick_params(colors='gray')
    ax3_pdf.spines['left'].set_color('gray')
    ax3_pdf.spines['bottom'].set_color('gray')
    ax3_pdf.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10,
    )
    l+=1  
    
    ax3_pdf.annotate(
        'Unprecedented\nlifetime\nexposure\nfor {} people'.format(str(int(np.round(ds_dmg['by_population_y0'].sel({'birth_year':2020,'lat':city_lat,'lon':city_lon},method='nearest').item(),-3)))),
        (1.1,0.6),
        xycoords=ax3_pdf.transAxes,
        fontsize=14,
        rotation='horizontal',
        color='gray',
        # fontweight='bold',
    )                   

    # City name
    ax3.annotate(
        '{}, {}'.format(city_name,cntry),
        (1960,ax3.get_ylim()[-1]),
        xycoords=ax3.transData,
        fontsize=16,
        rotation='horizontal',
        color='gray',
    )

    # axis labels ===================================================================

    # x axis label (time)
    x_i=1950
    y_i=-10
    x_f=2040
    y_f=y_i 
    con = ConnectionPatch(
        xyA=(x_i,y_i),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con)   

    con_arrow_top = ConnectionPatch(
        xyA=(x_f-2,y_f+1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_top)  

    con_arrow_bottom = ConnectionPatch(
        xyA=(x_f-2,y_f-1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_bottom) 
    ax.annotate(
        'Time',
        ((x_i+x_f)/2,y_f+1),
        xycoords=ax.transData,
        fontsize=12,
        color='gray',
    )

    # y axis label (Cumulative heatwave exposure since birth)
    x_i=1950
    y_i=-10
    x_f=x_i
    y_f=y_i + 61
    con = ConnectionPatch(
        xyA=(x_i,y_i),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con)   

    con_arrow_left = ConnectionPatch(
        xyA=(x_f-2,y_f-1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_left)  

    con_arrow_right = ConnectionPatch(
        xyA=(x_f+2,y_f-1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_right) 

    ax.annotate(
        'Cumulative heatwave exposure since birth',
        (x_i-10,(y_i+y_f)/5),
        xycoords=ax.transData,
        fontsize=12,
        rotation='vertical',
        color='gray',
    )

    # legend ===================================================================

    # bbox
    x0 = 1.5
    y0 = 0.5
    xlen = 0.5
    ylen = 0.5

    # space between entries
    legend_entrypad = 0.5

    # length per entry
    legend_entrylen = 0.75

    legend_font = 10
    legend_lw=2
    
    legendcols = list(colors.values())+['gray']+['lightgrey']
    handles = [
        Line2D([0],[0],linestyle='-',lw=legend_lw,color=legendcols[0]),
        Line2D([0],[0],linestyle='-',lw=legend_lw,color=legendcols[1]),
        Line2D([0],[0],linestyle='-',lw=legend_lw,color=legendcols[2]),
        Line2D([0],[0],linestyle='--',lw=legend_lw,color=legendcols[3]),
        Rectangle((0,0),1,1,color=legendcols[4]),
    ]
    labels= [
        '1.5 °C GMT warming by 2100',
        '2.5 °C GMT warming by 2100',
        '3.5 °C GMT warming by 2100',
        '99.99% pre-industrial \n lifetime exposure',
        'pre-industrial lifetime \n exposure histogram'
    ]
    ax.legend(
        handles, 
        labels, 
        bbox_to_anchor=(x0, y0, xlen, ylen), # bbox: (x, y, width, height)
        loc='upper left',
        ncol=1,
        fontsize=legend_font, 
        labelcolor='gray',
        mode="upper left", 
        borderaxespad=0.,
        frameon=False, 
        columnspacing=0.05, 
    )      

    # # population estimates
    # ds_dmg['population'].sel({'time':1990,'lat':city_lat,'lon':city_lon},method='nearest').sum(dim='age')

    # ds_dmg['by_population_y0'].sel({'birth_year':2020,'lat':city_lat,'lon':city_lon},method='nearest').item()
    
    # # getting estimate of all birth years that emerge in 1.5 and 3.5 pathways and how many these cohorts sum to
    # valid_bys=da_test_city.birth_year.where(da_test_city.loc[{'GMT':0}].max(dim='time')>da_pic_city_9999)
    # y1 = valid_bys.min(dim='birth_year')
    # y2 = valid_bys.max(dim='birth_year')
  
    # unprecedented=ds_dmg['by_population_y0'].sel(birth_year=np.arange(y1,y2+1),lat=city_lat,lon=city_lon,method='nearest').sum(dim='birth_year').round().item()    
    # print('{} thousand unprecedented born in {} and later under pathway {}'.format(unprecedented/10**3,y1,0))
    
    # valid_bys=da_test_city.birth_year.where(da_test_city.loc[{'GMT':20}].max(dim='time')>da_pic_city_9999)
    # y1 = valid_bys.min(dim='birth_year')
    # y2 = valid_bys.max(dim='birth_year')
    # unprecedented=ds_dmg['by_population_y0'].sel(birth_year=np.arange(y1,y2+1),lat=city_lat,lon=city_lon,method='nearest').sum(dim='birth_year').round().item()    
    # print('{} thousand unprecedented born in {} and later under pathway {}'.format(unprecedented/10**3,y1,20))        

    f.savefig(scripts_dir+'/figures/assessment/BiCC2_Norway/f1_concept_{}_{}.png'.format(flags['version'],cntry),dpi=1000,bbox_inches='tight')

#%%-----------------------------------------------------------------#
# Cohorts fractions for heatwaves for a specific country            #
# based on Fig.2 of Grant et al.(2025-)                             #
#-------------------------------------------------------------------#

def plot_Norway_BiCC2_fig2(
    df_GMT_strj,
    ds_pf_gs,
    da_gs_popdenom,
    gdf_country_borders,
    sims_per_step,
    flags,
    df_countries,
):

    # plot characteristics
    plot_var='unprec_99.99'
    x=12
    y=7
    markersize=10
    col_cbticlbl = 'gray'   # colorbar color of tick labels
    col_cbtic = 'gray'   # colorbar color of ticks
    col_cbedg = 'gray'   # colorbar color of edge
    cb_ticlen = 3.5   # colorbar length of ticks
    cb_ticwid = 0.4   # colorbar thickness of ticks
    cb_edgthic = 0   # colorbar thickness of edges between colors   
    by=2020
    cmap_whole = plt.cm.get_cmap('Reds')
    levels_cmap = np.arange(0,1.01,0.05)
    colors = [cmap_whole(i) for i in levels_cmap[:-1]]
    cmap_list_frac = mpl.colors.ListedColormap(colors,N=len(colors))
    ticks = np.arange(0,101,10)
    levels = np.arange(0,101,5)
    norm=mpl.colors.BoundaryNorm(levels,ncolors=len(levels)-1)
    l = 0 # letter indexing
    gmt_indices_152535 = [20,10,0]
    map_letters = {20:'g',10:'f',0:'e'}

    gmt_legend={
        GMT_indices_plot[0]:'1.5',
        GMT_indices_plot[1]:'2.5',
        GMT_indices_plot[2]:'3.5',
    }     

    f = plt.figure(figsize=(x,y))    
    gs0 = gridspec.GridSpec(10,3)

    # box plots
    # ax0 = f.add_subplot(gs0[0:4,0:2]) 
    ax0 = f.add_subplot(gs0[0:5,0:2]) 

    # pop totals
    # ax1 = f.add_subplot(gs0[4:,0:2],sharex=ax0)
    ax1 = f.add_subplot(gs0[5:,0:2],sharex=ax0)

    gs0.update(hspace=0)

    # maps
    gs01 = gridspec.GridSpecFromSubplotSpec(
        3,
        1, 
        subplot_spec=gs0[0:9,2:],
    )
    ax01 = f.add_subplot(gs01[0],projection=ccrs.Robinson())
    ax11 = f.add_subplot(gs01[1],projection=ccrs.Robinson())
    ax21 = f.add_subplot(gs01[2],projection=ccrs.Robinson()) 
    
    pos00 = ax21.get_position()
    cax00 = f.add_axes([
        pos00.x0-0.05,
        pos00.y0-0.1,
        pos00.width*1.95,
        pos00.height*0.2
    ])

    # get data
    df_list_gs = []
    extr='heatwavedarea'

    with open(data_dir+'{}/{}/isimip_metadata_{}_{}_{}.pkl'.format(flags['version'],extr,extr,flags['gmt'],flags['rm']), 'rb') as file:
        d_isimip_meta = pk.load(file)              
    with open(data_dir+'{}/{}/gridscale_aggregated_pop_frac_{}.pkl'.format(flags['version'],extr,extr), 'rb') as file:
        ds_pf_gs_plot = pk.load(file)

        ds_pf_gs_plot = ds_pf_gs_plot.sel(country='Norway')
    
    da_p_gs_plot = ds_pf_gs_plot[plot_var].loc[{
        'GMT':GMT_indices_plot,
        'birth_year':sample_birth_years,
    }]

    # print(type(da_p_gs_plot))
    # print(da_p_gs_plot)
    # print("---------------")
    # print(type(da_gs_popdenom))
    # print(da_gs_popdenom)


    sims_per_step = {}
    for step in GMT_labels:
        sims_per_step[step] = []
        for i in list(d_isimip_meta.keys()):
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                sims_per_step[step].append(i)  
                
    # --------------------------------------------------------------------
    # pf boxplot time series
    for step in GMT_indices_plot:
        da_pf_gs_plot_step = da_p_gs_plot.loc[{'run':sims_per_step[step],'GMT':step}].fillna(0) / da_gs_popdenom * 100
        df_pf_gs_plot_step = da_pf_gs_plot_step.to_dataframe(name='pf').reset_index()
        df_pf_gs_plot_step['GMT_label'] = df_pf_gs_plot_step['GMT'].map(gmt_legend)       
        df_pf_gs_plot_step['hazard'] = extr
        df_list_gs.append(df_pf_gs_plot_step)
    df_pf_gs_plot = pd.concat(df_list_gs)

    colors = dict(zip(list(gmt_legend.values()),['steelblue','darkgoldenrod','darkred']))
    p = sns.boxplot(
        data=df_pf_gs_plot[df_pf_gs_plot['hazard']==extr],
        x='birth_year',
        y='pf',
        hue='GMT_label',
        palette=colors,
        showcaps=False,
        showfliers=False,
        boxprops={
            'linewidth':0,
            'alpha':0.5
        },        
        ax=ax0,
    )
    p.legend_.remove()                  
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)      
    ax0.tick_params(colors='gray')
    ax0.set_ylim(0,100)
    ax0.spines['left'].set_color('gray')
    ax0.spines['bottom'].set_color('gray')      
    ax0.set_ylabel('$\mathregular{CF_{Heatwaves}}$ [%]',color='gray',fontsize=14)
    ax0.set_xlabel('Birth year',color='gray',fontsize=14)       
    ax0.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10
    )    
    l+=1 

    # bbox
    x0 = 0.025
    y0 = 0.7
    xlen = 0.2
    ylen = 0.3

    # space between entries
    legend_entrypad = 0.5

    # length per entry
    legend_entrylen = 0.75

    legend_font = 10
    legend_lw=3.5   

    legendcols = list(colors.values())
    handles = [
        Rectangle((0,0),1,1,color=legendcols[0]),\
        Rectangle((0,0),1,1,color=legendcols[1]),\
        Rectangle((0,0),1,1,color=legendcols[2])
    ]

    labels= [
        '1.5 °C GMT warming by 2100',
        '2.5 °C GMT warming by 2100',
        '3.5 °C GMT warming by 2100',    
    ]

    ax0.legend(
        handles, 
        labels, 
        bbox_to_anchor=(x0, y0, xlen, ylen), 
        loc = 'upper left',
        ncol=1,
        fontsize=legend_font, 
        labelcolor='gray',
        mode="expand", 
        borderaxespad=0.,\
        frameon=False, 
        columnspacing=0.05, 
        handlelength=legend_entrylen, 
        handletextpad=legend_entrypad
    )   

    # --------------------------------------------------------------------
    # populations

    ax1.spines['right'].set_visible(False)
    # ax1.spines['top'].set_visible(False)      
    ax1.tick_params(colors='gray')
    # ax1.set_ylim(0,100)
    ax1.spines['top'].set_color('gray')   
    ax1.spines['left'].set_color('gray')
    ax1.spines['bottom'].set_color('gray')    

    incomegroups = df_countries['incomegroup'].unique()
    income_countries = {}
    for category in incomegroups:
        income_countries[category] = list(df_countries.index[df_countries['incomegroup']==category])

    heights={}
    for category in incomegroups:
        heights[category] = da_gs_popdenom.loc[{
            'birth_year':np.arange(1960,2021,10),'country':income_countries[category]
        }].sum(dim='country').values / 10**6
    testdf = pd.DataFrame(heights)    
    testdf['birth_year'] = np.arange(1960,2021,10)
    testdf = testdf.set_index('birth_year')             
    testdf['total'] = testdf.sum(axis=1)
    p1 = testdf['total'].plot(
        kind='bar',
        # column='total',
        # stacked=True,
        color='gray',      
        ax=ax1,
        legend=False,
        rot=0.5
    ) 
    ax1.invert_yaxis()
    ax1.set_xlabel('Birth year',color='gray',fontsize=14)
    ax1.set_ylabel('Global cohort \n totals [in millions]',color='gray',fontsize=14)  
    ax1.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10
    )  
    # plot cohort sizes as text inside bars
    for i,by in enumerate(testdf.index.values):
        ax1.text(
            x=i,
            y=np.round(testdf['total'].loc[by]) - 10,
            s=str(int(np.round(testdf['total'].loc[by]))),
            horizontalalignment='center',
            verticalalignment='center',
            # transform=ax1.transData,
            fontsize=8,
            color='white',
        )
    
   
    
    l+=1
    # --------------------------------------------------------------------
    # maps of pop frac emergence for countries at 1, 2 and 3 deg pathways

    da_p_gs_plot = ds_pf_gs[plot_var].loc[{
        'GMT':gmt_indices_152535,
        'birth_year':by,
    }]
    df_list_gs = []
    for step in gmt_indices_152535:
        da_p_gs_plot_step = da_p_gs_plot.loc[{'run':sims_per_step[step],'GMT':step}].mean(dim='run')
        da_p_gs_plot_step = da_p_gs_plot_step / da_gs_popdenom.loc[{'birth_year':by}] * 100
        df_p_gs_plot_step = da_p_gs_plot_step.to_dataframe(name='pf').reset_index()
        df_p_gs_plot_step = df_p_gs_plot_step.assign(GMT_label = lambda x: np.round(df_GMT_strj.loc[2100,x['GMT']],1).values.astype('str'))
        df_list_gs.append(df_p_gs_plot_step)
    df_p_gs_plot = pd.concat(df_list_gs)
    df_p_gs_plot['pf'] = df_p_gs_plot['pf'].fillna(0)  
    gdf = cp(gdf_country_borders.reset_index())
    gdf_p = cp(gdf_country_borders.reset_index())
    robinson = ccrs.Robinson().proj4_init

    for ax,step in zip((ax01,ax11,ax21),gmt_indices_152535):
        gdf_p['pf']=df_p_gs_plot['pf'][df_p_gs_plot['GMT']==step].values
        ax.add_feature(feature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='powderblue', facecolor='powderblue'))
        gdf_p.to_crs(robinson).plot(
            ax=ax,
            column='pf',
            cmap=cmap_list_frac,
            norm=norm,
            cax=cax00,
            zorder=2,
            rasterized=True,
        )

        gdf.to_crs(robinson).plot(
            ax=ax,
            color='none', 
            edgecolor='black',
            linewidth=0.25,
            zorder=3,
        )
        
        gdf_robinson_bounds = gdf_p.to_crs(robinson).total_bounds # (minx,miny,maxx,maxy)
        
        # ax.set_global() # checking to see if this makes projections fully showing antarctica
        
        ax.set_title(
            '{} °C'.format(gmt_legend[step]),
            loc='center',
            fontweight='bold',
            fontsize=12,
            color='gray',       
        )
        
        ax.set_title(
            map_letters[step],
            loc='left',
            fontweight='bold',
            fontsize=10
        )    

        # triangles showing connections between GMT step box plot to map panels ---------------
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        if step == gmt_indices_152535[0]:
            x_h=1 
        elif step == gmt_indices_152535[1]:
            x_h=0.95                      
        elif step == gmt_indices_152535[-1]:
            x_h=0.9
        y_h= df_pf_gs_plot[(df_pf_gs_plot['birth_year']==by)&(df_pf_gs_plot['GMT']==step)]['pf'].median() / 100
        x_m=0
        y_m_i=0
        y_m_f=1.0
        con_low = ConnectionPatch(
            xyA=(x_h,y_h),
            xyB=(x_m,y_m_i),
            coordsA=ax0.transAxes,
            coordsB=ax.transAxes,
            color='lightgray',
            alpha=0.5,
            zorder=5,
        )
        con_hi = ConnectionPatch(
            xyA=(x_h,y_h),
            xyB=(x_m,y_m_f),
            coordsA=ax0.transAxes,
            coordsB=ax.transAxes,
            color='lightgray',
            alpha=0.5,
            zorder=5,
        )   
        con_vert = ConnectionPatch(
            xyA=(x_m,y_m_i),
            xyB=(x_m,y_m_f),
            coordsA=ax.transAxes,
            coordsB=ax.transAxes,
            color='lightgray',
            alpha=0.5,
            zorder=5,
        )
        
        line_low = con_low.get_path().vertices
        line_hi = con_hi.get_path().vertices
        
        ax0.add_artist(con_low)  
        ax0.add_artist(con_hi) 
        
        tri_coords = np.stack(
            (con_low.get_path().vertices[0],con_low.get_path().vertices[-1],con_hi.get_path().vertices[-1]),
            axis=0,
        )
        triangle = plt.Polygon(tri_coords,ec='lightgray',fc='lightgray',alpha=0.5,zorder=10,clip_on=False)    
        ax0.add_artist(triangle)           
        
        # triangles showing connections between GMT step box plot to map panels ---------------        
        
    cb = mpl.colorbar.ColorbarBase(
        ax=cax00, 
        cmap=cmap_list_frac,
        norm=norm,
        orientation='horizontal',
        spacing='uniform',
        ticks=ticks,
        drawedges=False,
    )

    cb.set_label(
        '$\mathregular{CF_{Heatwaves}}$ for 2020 birth cohort [%]',
        fontsize=14,
        color='gray'
    )
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(
        labelcolor=col_cbticlbl,
        color=col_cbtic,
        length=cb_ticlen,
        width=cb_ticwid,
        direction='out'
    )   
    cb.outline.set_edgecolor(col_cbedg)
    cb.outline.set_linewidth(cb_edgthic)   
    cax00.xaxis.set_label_position('top')   

    f.savefig(scripts_dir+'/figures/assessment/BiCC2_Norway/f2_combined_plot_popsizes.png',dpi=1000)


#%%-----------------------------------------------------------------#
# Conceptual plot for emergence in one location                     #
# based on Fig.1 of Grant et al.(2025-)                             #
#-------------------------------------------------------------------#
   

def plot_website_fig1(
    da_cohort_size,
    countries_mask,
    countries_regions,
    d_isimip_meta,
    flags,
    df_life_expectancy_5,
):
    # get data
    cntry=countries
    city_name='Madrid'
    # concept_bys = np.arange(1960,2021,30)
    concept_bys = np.arange(1960,2021,1)
    print("Conceptual plot for "+cntry)
    da_smple_cht = da_cohort_size.sel(country=cntry) # cohort absolute sizes in sample country
    da_smple_cht_prp = da_smple_cht / da_smple_cht.sum(dim='ages') # cohort relative sizes in sample country
    da_cntry = xr.DataArray(
        np.in1d(countries_mask,countries_regions.map_keys(cntry)).reshape(countries_mask.shape),
        dims=countries_mask.dims,
        coords=countries_mask.coords,
    )
    da_cntry = da_cntry.where(da_cntry,drop=True)
    # weights for latitude (probably won't use but will use population instead)
    lat_weights = np.cos(np.deg2rad(da_cntry.lat))
    lat_weights.name = "weights"   
    # city coords  
    city_lat = 40.4168 
    city_lon = 3.7038 

    ds_spatial = xr.Dataset(
        data_vars={
            'cumulative_exposure': (
                ['run','GMT','birth_year','time','lat','lon'],
                np.full(
                    (len(list(d_isimip_meta.keys())),
                    len(GMT_indices_plot),
                    len(concept_bys),
                    len(year_range),
                    len(da_cntry.lat.data),
                    len(da_cntry.lon.data)),
                    fill_value=np.nan,
                ),
            ),
        },
        coords={
            'lat': ('lat', da_cntry.lat.data),
            'lon': ('lon', da_cntry.lon.data),
            'birth_year': ('birth_year', concept_bys),
            'time': ('time', year_range),
            'run': ('run', np.arange(1,len(list(d_isimip_meta.keys()))+1)),
            'GMT': ('GMT', GMT_indices_plot)
        }
    )

    # load demography pickle
    with open(data_dir+'{}/gridscale/gridscale_dmg_{}.pkl'.format(flags['version'],cntry), 'rb') as f:
        ds_dmg = pk.load(f)                 

    # loop over simulations
    for i in list(d_isimip_meta.keys()): 

        print('Loading simulation {} of {}'.format(i,len(d_isimip_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA = pk.load(f)
            
        # mask to sample country and reduce spatial extent
        da_AFA = da_AFA.where(ds_dmg['country_extent']==1,drop=True)
        
        for step in GMT_indices_plot:
            
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                
                da_AFA_step = da_AFA.reindex(
                    {'time':da_AFA['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
                ).assign_coords({'time':year_range})                     
                                    
                # simple lifetime exposure sum
                da_le = xr.concat(
                    [(da_AFA_step.loc[{'time':np.arange(by,ds_dmg['death_year'].sel(birth_year=by).item()+1)}].cumsum(dim='time') +\
                    da_AFA_step.sel(time=ds_dmg['death_year'].sel(birth_year=by).item()) *\
                    (ds_dmg['life_expectancy'].sel(birth_year=by).item() - np.floor(ds_dmg['life_expectancy'].sel(birth_year=by)).item()))\
                    for by in concept_bys],
                    dim='birth_year',
                ).assign_coords({'birth_year':concept_bys})
                
                da_le = da_le.reindex({'time':year_range})
                
                ds_spatial['cumulative_exposure'].loc[{
                    'run':i,
                    'GMT':step,
                    'birth_year':concept_bys,
                    'time':year_range,
                    'lat':ds_dmg['country_extent'].lat.data,
                    'lon':ds_dmg['country_extent'].lon.data,
                }] = da_le.loc[{
                    'birth_year':concept_bys,
                    'time':year_range,
                    'lat':ds_dmg['country_extent'].lat.data,
                    'lon':ds_dmg['country_extent'].lon.data,
                }]

    # mean for the city           
    da_test_city = ds_spatial['cumulative_exposure'].sel({'lat':city_lat,'lon':city_lon},method='nearest').mean(dim='run')
    da_test_city = da_test_city.rolling(time=5,min_periods=5).mean()
    print("da_test_city ", da_test_city)

    # standard deviation for the city
    da_test_city_std = ds_spatial['cumulative_exposure'].sel({'lat':city_lat,'lon':city_lon},method='nearest').std(dim='run')
    da_test_city_std = da_test_city_std.rolling(time=5,min_periods=5).mean()

    # fill in 1st 4 years with 1s
    # first for mean
    for by in da_test_city.birth_year.data:
        for step in GMT_indices_plot:
            da_test_city.loc[{'birth_year':by,'GMT':step,'time':np.arange(by,by+5)}] = da_test_city.loc[{'birth_year':by,'GMT':step}].min(dim='time')
    # then for std        
    for by in da_test_city_std.birth_year.data:
        for step in GMT_indices_plot:
            da_test_city_std.loc[{'birth_year':by,'GMT':step,'time':np.arange(by,by+5)}] = da_test_city_std.loc[{'birth_year':by,'GMT':step}].min(dim='time')        
                
    # load PIC pickles
    with open(data_dir+'{}/{}/{}/gridscale_le_pic_{}_{}.pkl'.format(flags['version'],flags['extr'],cntry,flags['extr'],cntry), 'rb') as f:
        ds_pic = pk.load(f)   
    with open(data_dir+'{}/{}/{}/gridscale_pic_qntls_{}_{}.pkl'.format(flags['version'],flags['extr'],cntry,flags['extr'],cntry), 'rb') as f:
        ds_pic_qntl = pk.load(f)

    # plotting city lat/lon pixel doesn't give smooth kde
    print("ds_pic ", ds_pic)
    df_pic_city = ds_pic['lifetime_exposure'].sel({'lat':city_lat,'lon':city_lon},method='nearest').to_dataframe().drop(columns=['lat','lon',])
    print("df_pic_city ", df_pic_city)        
    da_pic_city_9999 = ds_pic_qntl['99.99'].sel({'lat':city_lat,'lon':city_lon},method='nearest')  
    print("da_pic_city_9999 ", da_pic_city_9999)

    # concept figure
    # ------------------------------------------------------------------   
    
    # plot building
    from mpl_toolkits.axes_grid1 import inset_locator as inset
    plt.rcParams['patch.linewidth'] = 0.1
    plt.rcParams['patch.edgecolor'] = 'k'
    colors = dict(zip(GMT_indices_plot,['steelblue','darkgoldenrod','darkred']))
    x=5
    y=1
    l = 0
    gmt_legend={
        GMT_indices_plot[0]:'1.5',
        GMT_indices_plot[1]:'2.5',
        GMT_indices_plot[2]:'3.5',
    }

    f,ax = plt.subplots(
        figsize=(x,y)
        )

    # ------------------------------------------------------------------   
    # 1960 time series
    

    for step in GMT_indices_plot:
        da_test_city.loc[{'birth_year':1960,'GMT':step}].plot.line(
            ax=ax,
            color=colors[step],
            linewidth=1,
        )
        # bold line for emergence
        da = da_test_city.loc[{'birth_year':1960,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        da.plot.line(
            ax=ax,
            color=colors[step],
            linewidth=3,
            zorder=4,
        )
             
    end_year=1960+np.floor(df_life_expectancy_5.loc[1960,cntry])
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    ax.set_xticks(np.arange(1960,2031,10))
    ax.set_xticklabels([1960,None,1980,None,2000,None,2020,None])
    ax.set_yticks([0,5])
    ax.set_yticklabels([None,5])     
    ax.annotate(
        'Born in 1960',
        (1965,ax.get_ylim()[-1]+2),
        xycoords=ax.transData,
        fontsize=10,
        fontweight='bold',
        rotation='horizontal',
        color='gray',
    )    
    ax.set_title(None)
    ax.annotate(
        letters[l],
        (1960,ax.get_ylim()[-1]+2),
        xycoords=ax.transData,
        fontsize=10,
        rotation='horizontal',
        color='k',
        fontweight='bold',
    )    
    l+=1      
      
    ax.set_xlim(
        1960,
        end_year,
    )
    ax.set_ylim(
        0,
        da_pic_city_9999+1,
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)    
    ax.tick_params(colors='gray')
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_color('gray')
    ax.hlines(
        y=da_pic_city_9999, 
        xmin=1960, 
        xmax=da_test_city.loc[{'birth_year':1960}].time.max()+10, 
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )

    # 1960 pdf
    ax_pdf_l = end_year+5
    ax_pdf_b = -2
    ax_pdf_w = 20
    ax_pdf_h = ax.get_ylim()[-1]+2
    ax_pdf = ax.inset_axes(
        bounds=(ax_pdf_l, ax_pdf_b, ax_pdf_w, ax_pdf_h),
        transform=ax.transData,
    )
    sns.histplot(
        data=df_pic_city.round(),
        y='lifetime_exposure',
        color='lightgrey',
        discrete = True,
        ax=ax_pdf
    )
    ax_pdf.hlines(
        y=da_pic_city_9999, 
        xmin=0, 
        xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    for step in GMT_indices_plot:
        ax_pdf.hlines(
            y=da_test_city.loc[{'birth_year':1960,'GMT':step}].max(), 
            xmin=0, 
            xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
            colors=colors[step], 
            linewidth=1, 
            linestyle='-', 
            label=gmt_legend[step], 
            zorder=2
        )
    ax_pdf.spines['right'].set_visible(False)
    ax_pdf.spines['top'].set_visible(False)      
    ax_pdf.set_ylabel(None)
    ax_pdf.set_xlabel(None)
    ax_pdf.set_ylim(-2,ax.get_ylim()[-1])
    ax_pdf.tick_params(colors='gray')
    ax_pdf.spines['left'].set_color('gray')
    ax_pdf.spines['bottom'].set_color('gray')
    ax_pdf.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10,
    )
    l+=1
        
    # ------------------------------------------------------------------       
    # 1990 time series
    ax2_l = 1990
    ax2_b = da_pic_city_9999 *2
    ax2_w = np.floor(df_life_expectancy_5.loc[1990,cntry])
    ax2_h = np.round(da_test_city.loc[{'birth_year':1990,'GMT':GMT_indices_plot[-1]}].max())
    ax2 = ax.inset_axes(
        bounds=(ax2_l, ax2_b, ax2_w, ax2_h),
        transform=ax.transData,
    )

    for step in GMT_indices_plot:
        da_test_city.loc[{'birth_year':1990,'GMT':step}].plot.line(
            ax=ax2,
            color=colors[step],
            linewidth=1,
        )
        # bold line for emergence
        da = da_test_city.loc[{'birth_year':1990,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        da.plot.line(
            ax=ax2,
            color=colors[step],
            linewidth=3,
            zorder=4,
        )    
                  
    end_year=1990+np.floor(df_life_expectancy_5.loc[1990,cntry])
    ax2.set_ylabel(None)
    ax2.set_xlabel(None)
    ax2.set_yticks([0,5,10])
    ax2.set_yticklabels([None,5,10])  
    ax2.set_xticks(np.arange(1990,2071,10))      
    ax2.set_xticklabels([None,2000,None,2020,None,2040,None,2060,None])
    ax2.set_xlim(
        1990,
        end_year,
    )
    ax2.set_ylim(
        0,
        np.round(da_test_city.loc[{'birth_year':1990,'GMT':GMT_indices_plot[-1]}].max())+1,
    )
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)  
    ax2.spines['left'].set_position(('data',1990)) 
    ax2.tick_params(colors='gray')
    ax2.spines['left'].set_color('gray')
    ax2.spines['bottom'].set_color('gray')
    
    ax2.annotate(
        'Born in 1990',
        (1995,ax2.get_ylim()[-1]),
        xycoords=ax2.transData,
        fontsize=10,
        fontweight='bold',
        rotation='horizontal',
        color='gray',
    )
    ax2.set_title(None)
    ax2.annotate(
        letters[l],
        (1990,ax2.get_ylim()[-1]),
        xycoords=ax2.transData,
        fontsize=10,
        rotation='horizontal',
        color='k',
        fontweight='bold',
    )     
    l+=1           

    # get time of first line to cross PIC thresh

    # does not work for some countries #

    emergences = []
    for step in GMT_indices_plot:
        da = da_test_city.loc[{'birth_year':1990,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        if np.any(da.notnull()): # testing for median
            da_t = da.time.where(da == da.min()).dropna(dim='time').item()
            emergences.append(da_t)
    print("da ", da)
    print("emeregences ", emergences)
    first_emerge = np.min(emergences)

    ax2.hlines(
        y=da_pic_city_9999, 
        xmin=first_emerge, 
        xmax=end_year, 
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )        

    # 1990 pdf

    ax2_pdf_l = end_year+5
    ax2_pdf_b = -2
    ax2_pdf_w = 20
    ax2_pdf_h = ax2.get_ylim()[-1]+2
    ax2_pdf = ax2.inset_axes(
        bounds=(ax2_pdf_l, ax2_pdf_b, ax2_pdf_w, ax2_pdf_h),
        transform=ax2.transData,
    )
    sns.histplot(
        data=df_pic_city.round(),
        y='lifetime_exposure',
        color='lightgrey',
        discrete = True,
        ax=ax2_pdf
    )
    ax2_pdf.hlines(
        y=da_pic_city_9999, 
        xmin=0, 
        xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    for step in GMT_indices_plot:
        ax2_pdf.hlines(
            y=da_test_city.loc[{'birth_year':1990,'GMT':step}].max(), 
            xmin=0, 
            xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
            colors=colors[step], 
            linewidth=1, 
            linestyle='-', 
            label=gmt_legend[step], 
            zorder=2
        )
    ax2_pdf.spines['right'].set_visible(False)
    ax2_pdf.spines['top'].set_visible(False)      
    ax2_pdf.set_ylabel(None)
    ax2_pdf.set_xlabel(None)
    ax2_pdf.set_ylim(-2,ax2.get_ylim()[-1])
    ax2_pdf.tick_params(colors='gray')
    ax2_pdf.spines['left'].set_color('gray')
    ax2_pdf.spines['bottom'].set_color('gray')
    ax2_pdf.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10,
    )
    l+=1   
    
    ax2_pdf.annotate(
        'Unprecedented\nlifetime\nexposure\nfor {} people'.format(str(int(np.round(ds_dmg['by_population_y0'].sel({'birth_year':1990,'lat':city_lat,'lon':city_lon},method='nearest').item(),-3)))),
        (1.1,0.3),
        xycoords=ax2_pdf.transAxes,
        fontsize=14,
        rotation='horizontal',
        color='gray',
        # fontweight='bold',
    )             

    # ------------------------------------------------------------------   
    # 2020 time series
    ax3_l = 2020
    ax3_b = np.round(da_test_city.loc[{'birth_year':1990,'GMT':GMT_indices_plot[-1]}].max()) * 1.5
    ax3_w = np.floor(df_life_expectancy_5.loc[2020,cntry])
    ax3_h = np.round(da_test_city.loc[{'birth_year':2020,'GMT':GMT_indices_plot[-1]}].max())
    ax3 = ax2.inset_axes(
        bounds=(ax3_l, ax3_b, ax3_w, ax3_h),
        transform=ax2.transData,
    )
    # plot mean lines
    for step in GMT_indices_plot:
        da_test_city.loc[{'birth_year':2020,'GMT':step}].plot.line(
            ax=ax3,
            color=colors[step],
            linewidth=1,
        )
        # bold line for emergence
        da = da_test_city.loc[{'birth_year':2020,'GMT':step}]
        da = da.where(da>da_pic_city_9999)
        da.plot.line(
            ax=ax3,
            color=colors[step],
            linewidth=3,
            zorder=4,
        )    

    end_year=2020+np.floor(df_life_expectancy_5.loc[2020,cntry])      
    
    ax3.set_ylabel(None)
    ax3.set_xlabel(None)
    ax3.set_yticks([0,5,10,15,20,25])
    ax3.set_yticklabels([None,5,10,15,20,25])   
    ax3.set_xticks(np.arange(2020,2101,10))      
    ax3.set_xticklabels([2020,None,2040,None,2060,None,2080,None,2100])
    ax3.set_xlim(
        2020,
        end_year,
    )
    ax3.set_ylim(
        0,
        np.round(da_test_city.loc[{'birth_year':2020,'GMT':GMT_indices_plot[-1]}].max())+1,
    )
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)  
    ax3.spines['left'].set_position(('data',2020))  
    ax3.tick_params(colors='gray')
    ax3.spines['left'].set_color('gray')
    ax3.spines['bottom'].set_color('gray')

    # get time of first line to cross PIC thresh

    # does not work for some countries 
    # emergences = []
    # for step in GMT_indices_plot:
    #     da = da_test_city.loc[{'birth_year':2020,'GMT':step}]
    #     da = da.where(da>da_pic_city_9999)
    #     if np.any(da.notnull()):
    #         da_t = da.time.where(da == da.min()).dropna(dim='time').item()
    #         emergences.append(da_t)
    # first_emerge = np.min(emergences)

    # ax3.hlines(
    #     y=da_pic_city_9999, 
    #     xmin=first_emerge, 
    #     xmax=end_year, 
    #     colors='grey', 
    #     linewidth=1, 
    #     linestyle='--', 
    #     label='99.99%', 
    #     zorder=1
    # )
    ax3.annotate(
        'Born in 2020',
        (2025,ax3.get_ylim()[-1]),
        xycoords=ax3.transData,
        fontsize=10,
        fontweight='bold',
        rotation='horizontal',
        color='gray',
    )
    ax3.set_title(None)
    ax3.annotate(
        letters[l],
        (2020,ax3.get_ylim()[-1]),
        xycoords=ax3.transData,
        fontsize=10,
        rotation='horizontal',
        color='k',
        fontweight='bold',
    ) 
    l+=1      

    # 2020 pdf
    ax3_pdf_l = end_year+5
    ax3_pdf_b = -2
    ax3_pdf_w = 20
    ax3_pdf_h = ax3.get_ylim()[-1]+2
    ax3_pdf = ax3.inset_axes(
        bounds=(ax3_pdf_l, ax3_pdf_b, ax3_pdf_w, ax3_pdf_h),
        transform=ax3.transData,
    )
    sns.histplot(
        data=df_pic_city.round(),
        y='lifetime_exposure',
        color='lightgrey',
        discrete = True,
        ax=ax3_pdf
    )
    ax3_pdf.hlines(
        y=da_pic_city_9999, 
        xmin=0, 
        xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
        colors='grey', 
        linewidth=1, 
        linestyle='--', 
        label='99.99%', 
        zorder=1
    )
    for step in GMT_indices_plot:
        ax3_pdf.hlines(
            y=da_test_city.loc[{'birth_year':2020,'GMT':step}].max(), 
            xmin=0, 
            xmax=df_pic_city['lifetime_exposure'][df_pic_city['lifetime_exposure']==0].count(),
            colors=colors[step], 
            linewidth=1, 
            linestyle='-', 
            label=gmt_legend[step], 
            zorder=2
        )
    ax3_pdf.spines['right'].set_visible(False)
    ax3_pdf.spines['top'].set_visible(False)      
    ax3_pdf.set_ylabel(None)
    ax3_pdf.set_xlabel(None)
    ax3_pdf.set_ylim(-2,ax3.get_ylim()[-1])
    ax3_pdf.tick_params(colors='gray')
    ax3_pdf.spines['left'].set_color('gray')
    ax3_pdf.spines['bottom'].set_color('gray')
    ax3_pdf.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10,
    )
    l+=1  
    
    # ax3_pdf.annotate(
    #     'Unprecedented\nlifetime\nexposure\nfor {} people'.format(str(int(np.round(ds_dmg['by_population_y0'].sel({'birth_year':2020,'lat':city_lat,'lon':city_lon},method='nearest').item(),-3)))),
    #     (1.1,0.6),
    #     xycoords=ax3_pdf.transAxes,
    #     fontsize=14,
    #     rotation='horizontal',
    #     color='gray',
    #     # fontweight='bold',
    # )                   

    # City name
    ax3.annotate(
        '{}, {}'.format(city_name,cntry),
        (1960,ax3.get_ylim()[-1]),
        xycoords=ax3.transData,
        fontsize=16,
        rotation='horizontal',
        color='gray',
    )

    # axis labels ===================================================================

    # x axis label (time)
    x_i=1950
    y_i=-10
    x_f=2040
    y_f=y_i 
    con = ConnectionPatch(
        xyA=(x_i,y_i),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con)   

    con_arrow_top = ConnectionPatch(
        xyA=(x_f-2,y_f+1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_top)  

    con_arrow_bottom = ConnectionPatch(
        xyA=(x_f-2,y_f-1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_bottom) 
    ax.annotate(
        'Time',
        ((x_i+x_f)/2,y_f+1),
        xycoords=ax.transData,
        fontsize=12,
        color='gray',
    )

    # y axis label (Cumulative heatwave exposure since birth)
    x_i=1950
    y_i=-10
    x_f=x_i
    y_f=y_i + 61
    con = ConnectionPatch(
        xyA=(x_i,y_i),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con)   

    con_arrow_left = ConnectionPatch(
        xyA=(x_f-2,y_f-1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_left)  

    con_arrow_right = ConnectionPatch(
        xyA=(x_f+2,y_f-1),
        xyB=(x_f,y_f),
        coordsA=ax.transData,
        coordsB=ax.transData,
        color='gray',
    )
    ax.add_artist(con_arrow_right) 

    ax.annotate(
        'Cumulative heatwave exposure since birth',
        (x_i-10,(y_i+y_f)/5),
        xycoords=ax.transData,
        fontsize=12,
        rotation='vertical',
        color='gray',
    )

    # legend ===================================================================

    # bbox
    x0 = 1.5
    y0 = 0.5
    xlen = 0.5
    ylen = 0.5

    # space between entries
    legend_entrypad = 0.5

    # length per entry
    legend_entrylen = 0.75

    legend_font = 10
    legend_lw=2
    
    legendcols = list(colors.values())+['gray']+['lightgrey']
    handles = [
        Line2D([0],[0],linestyle='-',lw=legend_lw,color=legendcols[0]),
        Line2D([0],[0],linestyle='-',lw=legend_lw,color=legendcols[1]),
        Line2D([0],[0],linestyle='-',lw=legend_lw,color=legendcols[2]),
        Line2D([0],[0],linestyle='--',lw=legend_lw,color=legendcols[3]),
        Rectangle((0,0),1,1,color=legendcols[4]),
    ]
    labels= [
        '1.5 °C GMT warming by 2100',
        '2.5 °C GMT warming by 2100',
        '3.5 °C GMT warming by 2100',
        '99.99% pre-industrial \n lifetime exposure',
        'pre-industrial lifetime \n exposure histogram'
    ]
    ax.legend(
        handles, 
        labels, 
        bbox_to_anchor=(x0, y0, xlen, ylen), # bbox: (x, y, width, height)
        loc='upper left',
        ncol=1,
        fontsize=legend_font, 
        labelcolor='gray',
        mode="upper left", 
        borderaxespad=0.,
        frameon=False, 
        columnspacing=0.05, 
    )      

    # # population estimates
    # ds_dmg['population'].sel({'time':1990,'lat':city_lat,'lon':city_lon},method='nearest').sum(dim='age')

    # ds_dmg['by_population_y0'].sel({'birth_year':2020,'lat':city_lat,'lon':city_lon},method='nearest').item()
    
    # # getting estimate of all birth years that emerge in 1.5 and 3.5 pathways and how many these cohorts sum to
    # valid_bys=da_test_city.birth_year.where(da_test_city.loc[{'GMT':0}].max(dim='time')>da_pic_city_9999)
    # y1 = valid_bys.min(dim='birth_year')
    # y2 = valid_bys.max(dim='birth_year')
  
    # unprecedented=ds_dmg['by_population_y0'].sel(birth_year=np.arange(y1,y2+1),lat=city_lat,lon=city_lon,method='nearest').sum(dim='birth_year').round().item()    
    # print('{} thousand unprecedented born in {} and later under pathway {}'.format(unprecedented/10**3,y1,0))
    
    # valid_bys=da_test_city.birth_year.where(da_test_city.loc[{'GMT':20}].max(dim='time')>da_pic_city_9999)
    # y1 = valid_bys.min(dim='birth_year')
    # y2 = valid_bys.max(dim='birth_year')
    # unprecedented=ds_dmg['by_population_y0'].sel(birth_year=np.arange(y1,y2+1),lat=city_lat,lon=city_lon,method='nearest').sum(dim='birth_year').round().item()    
    # print('{} thousand unprecedented born in {} and later under pathway {}'.format(unprecedented/10**3,y1,20))        


    f.savefig(scripts_dir+'/figures/assessment/website/f1_concept_{}_{}.png'.format(flags['version'],cntry),dpi=1000)


