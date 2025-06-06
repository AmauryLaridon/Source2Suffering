# ---------------------------------------------------------------
# Functions to compute emergence of exposure from noise
# ----------------------------------------------------------------

#               
#%%  ----------------------------------------------------------------
# IMPORT AND PATH 
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
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import interpolate
import cartopy.crs as ccrs
from settings import *
ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins = init()


# testing on finding countries with ar6 region  
    med = 'Mediterranean'
    gdf_ar6 = gpd.read_file('./data/shapefiles/IPCC-WGI-reference-regions-v4.shp')
    gdf_ar6 = gdf_ar6.loc[gdf_ar6['Name']==med]
    gdf_ar6 = gdf_ar6.loc[:,['Name','geometry']]
    gdf_ar6 = gdf_ar6.rename(columns={'Name':'name'})
    gdf_country = gdf_country_borders.loc[:,'geometry']
    gdf_country = gdf_country.reset_index()
    gdf_country.intersects(gdf_ar6['geometry'].iloc[0])
    gdf_med_countries = gdf_country.loc[gdf_country.intersects(gdf_ar6['geometry'].iloc[0])]
    lat = grid_area.lat.values
    lon = grid_area.lon.values    
    countries_med_3D = rm.mask_3D_geopandas(gdf_med_countries,lon,lat)
    ar6_regs_3D = rm.defined_regions.ar6.land.mask_3D(lon,lat)
    med_3D = ar6_regs_3D.isel(region=(ar6_regs_3D.names == 'Mediterranean')).squeeze()
    c_valid = []
    for c in gdf_med_countries.index:
        # next line gives nans outside country, 
        # 1 in parts of country in AR6 Medit
        # and 0 in parts outside Medit.
        c_in_med = med_3D.where(countries_med_3D.sel(region=c)==1) 
        c_area_in_med = c_in_med.weighted(grid_area/10**6).sum(dim=('lat','lon'))
        # c_area_out_med = xr.where(c_in_med==0,1,0).weighted(grid_area/10**6).sum(dim=('lat','lon'))
        c_area = countries_med_3D.sel(region=c).weighted(grid_area/10**6).sum(dim=('lat','lon'))
        c_area_frac = c_area_in_med.item() / c_area.item()
        if c_area_frac > 0.5:
            c_valid.append(c)
    countries_med_3D = countries_med_3D.loc[{'region':c_valid}]
    gdf_med_countries = gdf_med_countries.loc[c_valid]
    
    da_med_p = ds_pf_strj['unprec_country_b_y0'].loc[{'country':gdf_med_countries['name'].values}].mean(dim='run')
    da_med_pf = da_med_p.sum(dim='country') / ds_cohorts['by_population_y0'].loc[{'country':da_med_p.country.data}].sum(dim='country')

    
    ds_e_test = calc_exposure_trends_test(
        d_isimip_meta,
        grid_area,
        gdf_country_borders,
        flags,
    )
        
    lat = grid_area.lat.values
    lon = grid_area.lon.values

    # 3d mask for ar6 regions
    ar6_regs_3D = rm.defined_regions.ar6.land.mask_3D(lon,lat)

    # 3d mask for countries
    # countries_3D = rm.defined_regions.natural_earth_v5_0_0.countries_110.mask_3D(lon,lat) opting to use same geodataframe as analysis instead of regionmask
    countries_3D = rm.mask_3D_geopandas(gdf_country_borders.reset_index(),lon,lat)

    # temporary read in of absolute population emerging (y0 one, will have this in ds_pop_frac in future iterations)
    # ds_pf_strj = ds_pf_strj.assign_coords({'country':ds_cohorts.country.data})
    # ds_pf_strj['unprec_country_b_y0'] = xr.DataArray(
    #     data=np.full(
    #         (len(list(d_isimip_meta.keys())),len(ds_cohorts.country.data),len(birth_years),len(GMT_labels)),
    #         fill_value=np.nan
    #     ),
    #     dims=['run','country','birth_year','GMT'],
    #     coords={
    #             'run': ('run', list(d_isimip_meta.keys())),
    #             'birth_year': ('birth_year', birth_years),
    #             'country': ('country', ds_cohorts.country.data),
    #             'GMT': ('GMT', GMT_labels),
    #         }
    # )
    # dataset for emergence masks for summing across runs
    ds_testing = xr.Dataset(
        data_vars={
            'emergence_masks': (
                ['run','country','birth_year','GMT'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(ds_cohorts.country.data),len(birth_years),len(GMT_labels)),
                    fill_value=np.nan,
                ),
            ),
            'totals_emergence_masks': (
                ['country','birth_year','GMT'],
                np.full(
                    (len(ds_cohorts.country.data),len(birth_years),len(GMT_labels)),
                    fill_value=np.nan,
                ),
            ),
            'unprec_country_b_y0': (
                ['run','country','birth_year','GMT'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(ds_cohorts.country.data),len(birth_years),len(GMT_labels)),
                    fill_value=np.nan,
                ),
            ),       
            'mean_unprec_country_b_y0': (
                ['country','birth_year','GMT'],
                np.full(
                    (len(ds_cohorts.country.data),len(birth_years),len(GMT_labels)),
                    fill_value=np.nan,
                ),
            ),       
        },
        coords={
            'run': ('run', list(d_isimip_meta.keys())),
            'birth_year': ('birth_year', birth_years),
            'country': ('country', ds_cohorts.country.data),
            'GMT': ('GMT', GMT_labels),
            'time': ('time', year_range),
        }    
    )

    for i in list(d_isimip_meta.keys()):
        s=0
        for step in GMT_labels:
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                # load pickle
                with open('./data/pickles/da_emergence_mask_{}_{}_{}_{}_{}.pkl'.format(flags['gmt'],flags['extr'],flags['rm'],i,step), 'rb') as f:
                    da_emergence_mask_run_step = pk.load(f)
                da_birthyear_emergence_mask = xr.where(da_emergence_mask_run_step.sum(dim='time')>0,1,0)
                ds_testing['emergence_masks'].loc[{
                    'country':ds_cohorts.country.data,
                    'birth_year':birth_years,
                    'run':i,
                    'GMT':step,
                }] = da_birthyear_emergence_mask
                ds_testing['unprec_country_b_y0'].loc[{
                    'country':ds_cohorts.country.data,
                    'birth_year':birth_years,
                    'run':i,
                    'GMT':step,            
                }] = ds_cohorts['by_population_y0'].where(da_birthyear_emergence_mask==1)             
                s+=1
                
    ds_testing['totals_emergence_masks'] = ds_testing['emergence_masks'].sum(dim='run')

    # testing 2020, step 28 bar plots of pf and ae per country, and then aggregate pf and ae for each sim
    step=28
    by=2020
    sims = []
    for i in list(d_isimip_meta.keys()):
        if d_isimip_meta[i]['GMT_strj_valid'][step]:
            sims.append(i)

    pfs = [] # sample for this GMT, by to check mean against heatmap
    for n,i in enumerate(list(d_isimip_meta.keys())):
        if d_isimip_meta[i]['GMT_strj_valid'][step]:
            
            # initiate plotting axes
            f,(ax1,ax2,ax3,ax4) = plt.subplots(
                nrows=4,
                ncols=1,
                figsize=(10,10),
            )
            
            # emergence mask
            em = ds_testing['emergence_masks'].loc[{'birth_year':by,'run':i,'GMT':step}].where(ds_testing['emergence_masks'].loc[{'birth_year':by,'run':i,'GMT':step}]!=0).notnull().assign_coords({'country':range(len(ds_cohorts.country.data))})
            
            # x axis labels of country strings based on em
            x_tick_labels = ds_cohorts.country.assign_coords({'country':range(len(ds_ae_strj.country.data))}).where(em,drop=True)
            
            # plot ae per country (ds_ae_strj)
            p_ae = ds_ae_strj['age_emergence'].loc[{'run':i,'GMT':step,'birth_year':by}].assign_coords({'country':range(len(ds_cohorts.country.data))})
            p_ae = p_ae.where(em,drop=True)
            weights = ds_cohorts['by_y0_weights'].loc[{'birth_year':by}].assign_coords({'country':range(len(ds_ae_strj.country.data))}).where(em,drop=True)
            mean_ae = p_ae.weighted(weights).mean(dim='country').item()
            ax1.bar(x_tick_labels,p_ae.values)
            ax1.set_title(
                'mean ae: {}'.format(str(int(np.round(mean_ae,0)))),
                loc='center',
                fontweight='bold',
            )     
            ax1.set_title(
                '{} \n{} \n{} \n{}°C @ 2100 \n{} birth cohort'.format(d_isimip_meta[i]['model'].split('/')[-1],d_isimip_meta[i]['gcm'],d_isimip_meta[i]['rcp'],np.round(df_GMT_strj.loc[2100,step],1),by),
                loc='right',
                fontweight='bold'
            )
            ax1.set_ylabel(
                'age emergence', 
                va='center', 
                rotation='vertical',
                labelpad=10,
            )            
            
            # plot p per country  ds_testing['unprec_per_country_b_y0']
            p_p = ds_testing['unprec_country_b_y0'].loc[{'run':i,'GMT':step,'birth_year':by}].assign_coords({'country':range(len(ds_cohorts.country.data))}).where(em,drop=True) * 1000
            sum_p = p_p.sum(dim='country')
            ax2.bar(x_tick_labels,p_p.values)
            ax2.set_title(
                'sum p: {} million'.format(str(int(np.round(sum_p,0) / 10**6))),
                loc='center',
                fontweight='bold',
            )
            ax2.set_ylabel(
                'population emergence', 
                va='center', 
                rotation='vertical',
                labelpad=10,
            )
            
            # plot pf per country (each country is frac of global pop; ensure closes to 1) ds_pf_strj
            p_pf = p_p / (ds_cohorts['by_population_y0'].loc[{'birth_year':by}].assign_coords({'country':range(len(ds_cohorts.country.data))}).sum(dim='country') * 1000)
            p_pf = p_pf.where(em,drop=True)
            pf = ds_pf_strj['frac_unprec_all_b_y0'].loc[{'run':i,'GMT':step,'birth_year':by}].item()
            pfs.append(pf)
            ax3.bar(x_tick_labels,p_pf.values)
            ax3.set_title(
                'pf: {}'.format(str(np.round(pf,2))),
                loc='center',
                fontweight='bold',
            )
            ax3.set_ylabel(
                'population fraction \n emerged', 
                va='center', 
                rotation='vertical',
                labelpad=10,
            )
            
            # plot lifetime emergence and pic threshold
            p_le = ds_le['lifetime_exposure'].loc[{'run':i,'GMT':step,'birth_year':by}].assign_coords({'country':range(len(ds_cohorts.country.data))})
            mean_le = p_le.weighted(weights).mean(dim='country').item()
            p_le = p_le.where(em,drop=True)
            p_pic = ds_exposure_pic['mmm_pic'].assign_coords({'country':range(len(ds_cohorts.country.data))}).where(em,drop=True)
            ax4.bar(x_tick_labels,p_le.values)
            ax4.plot(
                x_tick_labels,
                p_pic.values,
                marker='o',
                linestyle='',
                color='r'
            )
            ax4.set_title(
                'mean le: {}'.format(str(np.round(mean_le,2))),
                loc='center',
                fontweight='bold',
            )
            ax4.set_ylabel(
                'lifetime exposure', 
                va='center', 
                rotation='vertical',
                labelpad=10,
            )
            ax4.set_xticklabels(x_tick_labels.values, rotation='vertical')
            
            # ax stuff
            for n,ax in enumerate((ax1,ax2,ax3,ax4)):
                ax.set_title(
                    letters[n],
                    loc='left',
                    fontweight='bold',
                )
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)                 
                if n < 3:
                    ax.tick_params(labelbottom=False)
                    
            f.savefig('./figures/testing/ae_p_pf_le_GMT_{}_{}_{}_{}_{}'.format(step,d_isimip_meta[i]['model'].split('/')[-1],d_isimip_meta[i]['gcm'],d_isimip_meta[i]['rcp'],flags['extr']),dpi=800)

    # mean pf
    mean_pf = np.mean(pfs)

    # scatter plots of country-mean ae and global pf per run
    gmts2100 = np.round(df_GMT_strj.loc[2100,[0,5,10,15,20,25]].values,1)

    for by in sample_birth_years:
        for step in GMT_indices:
            
            # initiate plotting axes
            f,((ax1,ax2),(ax3,ax4)) = plt.subplots(
                nrows=2,
                ncols=2,
                figsize=(10,7),
            )

            # age emergence in ax 1 and 2
            ds_plt = ds_ae_strj['age_emergence']                             
            ds_plt_gmt = ds_plt.loc[{'birth_year':by}]
            ds_plt_gmt = ds_plt_gmt.weighted(ds_cohorts['by_y0_weights'].loc[{'birth_year':by}]).mean(dim='country')
            p = ds_plt_gmt.to_dataframe().reset_index(level="run")
            x = p.index.values
            y = p['age_emergence'].values
            ax1.scatter(
                x,
                y,
            )
            ax1.plot(
                GMT_labels,
                ds_plt_gmt.mean(dim='run').values,
                marker='o',
                linestyle='',
                color='r'
            )
            ax1.set_title(
                '{} birth cohort'.format(str(by)),
                loc='center',
                fontweight='bold',
            )
            ax1.set_ylabel(
                'age emergence', 
                va='center', 
                rotation='vertical',
                labelpad=10,
            )                                               
            ax1.set_xticks(
                ticks=[0,5,10,15,20,25],
                labels=None,
            )

            ds_plt_by = ds_plt.loc[{'GMT':step}]
            ds_plt_by = ds_plt_by.weighted(ds_cohorts['by_y0_weights']).mean(dim='country')
            p = ds_plt_by.to_dataframe().reset_index(level="run")
            x = p.index.values
            y = p['age_emergence'].values
            ax2.scatter(
                x,
                y,
            )
            ax2.plot(
                birth_years,
                ds_plt_by.mean(dim='run').values,
                marker='o',
                linestyle='',
                color='r'
            )   
            ax2.set_title(
                '{} @ 2100 [°C]'.format(str(np.round(df_GMT_strj.loc[2100,step],1))),
                loc='center',
                fontweight='bold',
            )       

            # pf in ax 3 and 4
            ds_plt = ds_pf_strj['frac_unprec_all_b_y0']                                                   
            ds_plt_gmt = ds_plt.loc[{'birth_year':by}]
            p = ds_plt_gmt.to_dataframe().reset_index(level="run")
            x = p.index.values
            y = p['frac_unprec_all_b_y0'].values
            ax3.scatter(
                x,
                y,
            )
            ax3.plot(
                GMT_labels,
                ds_plt_gmt.mean(dim='run').values,
                marker='o',
                linestyle='',
                color='r'
            )
            ax3.set_ylabel(
                'population fraction', 
                va='center', 
                rotation='vertical',
                labelpad=10,
            )          
            ax3.set_xlabel(
                'GMT anomaly at 2100 [°C]', 
                va='center', 
                labelpad=10,
            )                                           
            ax3.set_xticks(
                ticks=[0,5,10,15,20,25],
                labels=gmts2100,
            )

            ds_plt_by = ds_plt.loc[{'GMT':step}]
            p = ds_plt_by.to_dataframe().reset_index(level="run")
            x = p.index.values
            y = p['frac_unprec_all_b_y0'].values
            ax4.scatter(
                x,
                y,
            )
            ax4.plot(
                birth_years,
                ds_plt_by.mean(dim='run').values,
                marker='o',
                linestyle='',
                color='r'
            )   
            ax4.set_xlabel(
                'Birth year', 
                va='center', 
                labelpad=10,
            )         
                    
            # ax stuff
            for n,ax in enumerate((ax1,ax2,ax3,ax4)):
                ax.set_title(
                    letters[n],
                    loc='left',
                    fontweight='bold',
                )
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)                 
                if n < 3:
                    ax.tick_params(labelbottom=False)        
            
            f.savefig('./figures/testing/ae_pf_GMT_by_scatterplots_{}_{}_{}.png'.format(step,by,flags['extr']),dpi=800)
            
            plt.show()
            
            
            
    # heatmaps of population totals
    # ['run','country','birth_year','GMT'],

    ds_testing['unprec_total'] = ds_testing['unprec_country_b_y0'].sum(dim='country')
    ds_testing['unprec_total'] = ds_testing['unprec_total'].where(ds_testing['unprec_total']!=0) * 1000 / 10**6 # with unit conversion
    ds_testing['unprec_mean'] = ds_testing['unprec_total'].mean(dim='run')
    gmts2100 = np.round(df_GMT_strj.loc[2100,[0,5,10,15,20,25]].values,1)

    ds_testing['unprec_mean'].loc[{
        'birth_year':np.arange(1960,2021)
    }].plot(
        x='birth_year',
        y='GMT',
        add_colorbar=True,
        levels=10,
        cbar_kwargs={
            'label':'Population totals'
        }
    ) 
    p.axes.set_yticks(
        ticks=[0,5,10,15,20,25],
        labels=gmts2100
    )
    p.axes.set_xticks(
        ticks=np.arange(1960,2025,10),
    )    
    p.axes.set_ylabel('GMT anomaly at 2100 [°C]')
    p.axes.set_xlabel('Birth year')
    p.axes.figure.savefig('./figures/pop_by_heatmap_{}_{}_{}.png'.format(flags['extr'],flags['gmt'],flags['rm']))    
    plt.show()



            # # axis stuff
            # for n,ax in enumerate((ax1,ax2,ax3,ax4)):
            #     ax.spines['right'].set_visible(False)
            #     ax.spines['top'].set_visible(False)            
    # get fraction of emergence totals (i.e. per GMT and birth year, what are the fraction of sims that have emergence occuring)
    # need sum divided by number of non-nan runs per GMT, birth year
    test = xr.where(ds_testing['emergence_masks'].sum(dim=['run','birth_year','country'])>0,1,0)
    ds_testing = ds_testing.assign_coords({'country':range(len(ds_cohorts.country.data))})


    for testyear in np.arange(year_start,year_ref+1,20):
        
        # checking age emergence cohort
        test1 = ds_ae_strj['age_emergence'].loc[{'birth_year':testyear}].mean(dim='run')
        test1 = test1.assign_coords({'country':range(len(ds_ae_strj.country.data))})
        p1 = test1.plot(figsize=(12,12))
        p1.axes.figure.savefig('./figures/testing/p1_ae_{}_heatmap_{}.png'.format(testyear,flags['extr']),dpi=500)
        plt.show()

        # compare p1 with same dimension heatmap for flood trends
        test4_trend = ds_e['mean_exposure_trend_country_time_ranges'].loc[{'year':testyear}].assign_coords({'country':range(len(ds_ae_strj.country.data))})
        p4 = test4_trend.where(test4_trend!=0).plot(x='GMT',y='country',figsize=(12,12),cmap='RdBu',levels=20)
        p4.axes.figure.savefig('./figures/testing/p4_trends_{}_heatmap_{}.png'.format(testyear,flags['extr']),dpi=500)
        plt.show()
        
        p5 = ds_testing['totals_emergence_masks'].loc[{'birth_year':testyear}].plot(figsize=(12,12))
        p5.axes.figure.savefig('./figures/testing/p5_emergencetotals_{}_heatmap_{}.png'.format(testyear,flags['extr']),dpi=500)
        plt.show()
        
        p6 = ds_testing['mean_unprec_country_b_y0'].loc[{'birth_year':testyear}].plot(figsize=(12,12))
        p6.axes.figure.savefig('./figures/testing/p6_unprecmean_{}_heatmap_{}.png'.format(testyear,flags['extr']),dpi=500)
        plt.show()
        
            
    # # checking valid runs per GMT level:
    sim_counts = []
    for step in GMT_labels:
        print('step {}'.format(step))
        c=0
        for i in list(d_isimip_meta.keys()):
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                c+=1
        sim_counts.append(c)
        # print('step {} has {} runs included'.format(step,c))
    da_sim_counts = xr.DataArray(
        data=sim_counts,
        coords={'GMT':ds_pf_strj.GMT.data},
    )
    p_sc = da_sim_counts.plot(marker='o')
    p_sc[0].axes.figure.savefig('./figures/testing/p_sc_{}'.format(flags['extr']),dpi=500)
        

        
    # checking year GMT mapping for last GMT step    
    step = 5
    s = 0
    for i in list(d_isimip_meta.keys()):
        if d_isimip_meta[i]['GMT_strj_valid'][step]:
            # load AFA data of that run
            with open('./data/pickles/isimip_AFA_{}_{}.pkl'.format(flags['extr'],str(i)), 'rb') as f:
                da_AFA = pk.load(f)          
            # da_AFA = da_AFA.reindex(
            #             {'time':da_AFA['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
            #         ).assign_coords({'time':year_range}) 
            da_AFA = da_AFA.reindex(
                        {'time':da_AFA['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
                    )
            plt.plot(year_range,da_AFA.time.data)
            s+=1
            plt.show()
    print(s)    

    # testing age emergence results with different order of ensemble aggregation (highly sensitive here)
    test1 = ds_ae_strj['age_emergence'].mean(dim=('run')) # first mean per run, then across countries (unfair weight to countries that infrequently emerge)
    test1 = test1.weighted(ds_cohorts['by_y0_weights']).mean(dim='country')
    test1.plot(x='birth_year',y='GMT',)   

    # first mean across countries, then runs (preferred method)
    test2 = ds_ae_strj['age_emergence'].weighted(ds_cohorts['by_y0_weights']).mean(dim='country')
    test2 = test2.mean(dim=('run')) # first mean per run, then across countries (unfair weight to countries that infrequently emerge)
    test2.plot(x='birth_year',y='GMT',)   

    # original (basically same as test2)
    test3 = ds_ae_strj['age_emergence'].weighted(ds_cohorts['by_y0_weights']).mean(dim=('country','run'))
    test3.plot(x='birth_year',y='GMT',)   

    for i in list(d_isimip_meta.keys()):
        with open('./data/pickles/isimip_AFA_{}_{}.pkl'.format(flags['extr'],str(i)), 'rb') as f:
                da_AFA = pk.load(f)
        for step in GMT_labels:
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                da_AFA = da_AFA.reindex(
                    {'time':da_AFA['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
                ).assign_coords({'time':year_range})
                da_AFA_country_weighted_sum = da_AFA.weighted(countries_3D*grid_area/10**6).sum(dim=('lat','lon'))
                y=1960
                da_AFA_country_weighted_sum = da_AFA_country_weighted_sum.loc[{'time':np.arange(y,y+81)}]
                stats_y_country = vectorize_lreg(da_AFA_country_weighted_sum)
                slope_y_country = stats_y_country[0]
                ds_e['exposure_trend_country_time_ranges'].loc[{
                    'run':i,
                    'GMT':step,
                    'country':countries_3D.region.data,
                    'year':y,
                }] = slope_y_country
                mean_floods = ds_e['exposure_trend_country_time_ranges'].loc[{
                    'run':i,
                    'GMT':step,
                    'country':countries_3D.region.data,
                    'year':y,
                }].mean(dim='country').item()
                print('simulation {}, GMT {} has country mean flood trends at {}'.format(i,step,mean_floods))        
                
    # getting max country
    max_list = slope_y_country.region.where(slope_y_country == slope_y_country.max(dim='region')).values
    max = max_list[~np.isnan(max_list)].item()
    ds_cohorts.country.data[int(max)]
        
#%% ----------------------------------------------------------------
# grid scale
# ------------------------------------------------------------------

# from gridscale import *

# if flags['gridscale']:
    
#     ds_le_gs, ds_ae_gs, ds_pf_gs = grid_scale_emergence(
#         d_isimip_meta,
#         d_pic_meta,
#         flags['extr'],
#         da_cohort_size,
#         countries_regions,
#         countries_mask,
#         df_life_expectancy_5,
#         GMT_indices,
#         da_population,
#     )
    
# else:
    
#     # load pickled aggregated lifetime exposure, age emergence and pop frac datasets
#     with open('./data/pickles/gridscale_aggregated_lifetime_exposure_{}.pkl'.format(flags['extr']), 'rb') as f:
#         ds_le_gs = pk.load(f)
#     with open('./data/pickles/gridscale_aggregated_age_emergence_{}.pkl'.format(flags['extr']), 'rb') as f:
#         ds_ae_gs = pk.load(f)
#     with open('./data/pickles/gridscale_aggregated_pop_frac_{}.pkl'.format(flags['extr']), 'rb') as f:
#         ds_pf_gs = pk.load(f)

# # load spatially explicit datasets
# d_gs_spatial = {}
# for cntry in sample_countries:
#     with open('./data/pickles/gridscale_spatially_explicit_{}_{}.pkl'.format(flags['extr'],cntry), 'rb') as f:
#         d_gs_spatial[cntry] = pk.load(f)

#----------------------------------------------------------------
# plot emergence stuff
# ------------------------------------------------------------------

# plot pop frac and age emergence across GMT for stylized trajectories
# top panel; (y: frac unprecedented, x: GMT anomaly @ 2100)
# bottom panel; (y: age emergence, x: GMT anomaly @ 2100)
# plots both approaches to frac unprecedented; exposed vs full cohorts
# issue in these plots that early birth years for high warming trajectories won't live till 3-4 degree warming, but we misleadingly plot along these points
    # so, for e.g. 1970 BY, we need to limit the line up to life expectancy
    # will need cohort weighted mean of life expectancy across countries
    # 
    
    # # add emergence mask since i forgot to do it in gridscale.py (added it there but haven't rerun 10 Jan)
    # for cntry in sample_countries:
    #     d_gs_spatial[cntry]['emergence_mask'] = (['run','GMT','birth_year','lat','lon'],np.full(
    #                     (len(list(d_isimip_meta.keys())),len(GMT_indices),len(sample_birth_years),len(d_gs_spatial[cntry].lat.data),len(d_gs_spatial[cntry].lon.data)),
    #                     fill_value=np.nan,
    #                 ))
    #     for i in list(d_isimip_meta.keys()):
    #         for step in GMT_indices:
    #             if os.path.isfile('./data/pickles/gridscale_exposure_mask_{}_{}_{}_{}.pkl'.format(flags['extr'],cntry,i,step)):
    #                 with open('./data/pickles/gridscale_exposure_mask_{}_{}_{}_{}.pkl'.format(flags['extr'],cntry,i,step), 'rb') as f:
    #                     da_birthyear_exposure_mask = pk.load(f)
    #                     d_gs_spatial[cntry]['emergence_mask'].loc[{'run':i,'GMT':step}] = da_birthyear_exposure_mask.loc[{'birth_year':sample_birth_years}]
    
    # cntry = 'Canada'
    # ind_cntry = countries_regions.map_keys(cntry)
    # mask = xr.DataArray(
    #     np.in1d(countries_mask,ind_cntry).reshape(countries_mask.shape),
    #     dims=countries_mask.dims,
    #     coords=countries_mask.coords,
    # )
        
    # for cntry in sample_countries:
    #     ind_cntry = countries_regions.map_keys(cntry)
    #     mask = xr.DataArray(
    #         np.in1d(countries_mask,ind_cntry).reshape(countries_mask.shape),
    #         dims=countries_mask.dims,
    #         coords=countries_mask.coords,
    #     )        
    #     for analysis in ['lifetime_exposure','age_emergence','emergence_mask','population_emergence']:   
    #         if analysis == 'emergence_mask':
    #             d_gs_spatial[cntry][analysis] = d_gs_spatial[cntry][analysis].where(mask)
    #         if cntry != 'Russian Federation':
    #             projection = ccrs.PlateCarree()
    #         else:
    #             projection = ccrs.LambertConformal(central_longitude=36.6, central_latitude=53.7, cutoff=30)
    #         p = d_gs_spatial[cntry][analysis].loc[{
    #             'birth_year':birth_years_plot,
    #             'GMT':GMT_indices_plot,
    #         }].mean(dim='run').plot(
    #             col='birth_year',
    #             row='GMT',
    #             transform=ccrs.PlateCarree(),
    #             subplot_kws={"projection": projection},
    #             aspect=2,
    #             size=3
    #         )
    #         for ax in p.axes.flat:
    #             ax.coastlines()
    #             ax.gridlines()
    #         p.fig.savefig('./figures/gridscale_sample_{}_{}.png'.format(cntry,analysis))
    
                        
    # # spatial lifetime exposure dataset (subsetting birth years and GMT steps to reduce data load) per country
    #     # can also add spatial age emergence to here
    # ds_gs_spatial = xr.Dataset(
    #     data_vars={
    #         'lifetime_exposure': (
    #             ['run','GMT','birth_year','lat','lon'],
    #             np.full(
    #                 (len(list(d_isimip_meta.keys())),len(GMT_indices_plot),len(birth_years_plot),len(countries_mask.lat.data),len(countries_mask.lon.data)),
    #                 fill_value=np.nan,
    #             ),
    #         ),
    #         'age_emergence': (
    #             ['run','GMT','birth_year','lat','lon'],
    #             np.full(
    #                 (len(list(d_isimip_meta.keys())),len(GMT_indices_plot),len(birth_years_plot),len(countries_mask.lat.data),len(countries_mask.lon.data)),
    #                 fill_value=np.nan,
    #             ),
    #         ),
    #         'population_emergence': (
    #             ['run','GMT','birth_year','lat','lon'],
    #             np.full(
    #                 (len(list(d_isimip_meta.keys())),len(GMT_indices_plot),len(birth_years_plot),len(countries_mask.lat.data),len(countries_mask.lon.data)),
    #                 fill_value=np.nan,
    #             ),
    #         ),
    #         'emergence_mask': (
    #             ['run','GMT','birth_year','lat','lon'],
    #             np.full(
    #                 (len(list(d_isimip_meta.keys())),len(GMT_indices_plot),len(birth_years_plot),len(countries_mask.lat.data),len(countries_mask.lon.data)),
    #                 fill_value=np.nan,
    #             ), 
    #         )
    #     },
    #     coords={
    #         'lat': ('lat', countries_mask.lat.data),
    #         'lon': ('lon', countries_mask.lon.data),
    #         'birth_year': ('birth_year', birth_years_plot),
    #         'run': ('run', np.arange(1,len(list(d_isimip_meta.keys()))+1)),
    #         'GMT': ('GMT', GMT_indices_plot)
    #     }
    # )
    
    # for cntry in sample_countries:
    #     for analysis in ['lifetime_exposure','age_emergence','emergence_mask']:
    #         ds_gs_spatial[analysis].loc[{
    #             'lat':d_gs_spatial[cntry].lat.data,
    #             'lon':d_gs_spatial[cntry].lon.data,
    #             'birth_year':birth_years_plot,
    #             'GMT':GMT_indices_plot,            
    #         }] = d_gs_spatial[cntry][analysis].loc[{
    #             'lat':d_gs_spatial[cntry].lat.data,
    #             'lon':d_gs_spatial[cntry].lon.data,
    #             'birth_year':birth_years_plot,
    #             'GMT':GMT_indices_plot,
    #         }]
    
    # ind_cntrs = []
    # for cntry in sample_countries:
    #     ind_cntrs.append(countries_regions.map_keys(cntry))
    # mask = xr.DataArray(
    #     np.in1d(countries_mask,ind_cntrs).reshape(countries_mask.shape),
    #     dims=countries_mask.dims,
    #     coords=countries_mask.coords,
    # )
    # ds_gs_spatial = ds_gs_spatial.where(mask)
    # plottable = ds_gs_spatial['lifetime_exposure'].mean(dim='run')
    # plottable.plot(transform=ccrs.PlateCarree(),col='birth_year',row='GMT',subplot_kws={'projection':ccrs.PlateCarree()})
    
    
    # for analysis in ['lifetime_exposure','age_emergence','emergence_mask']:
        
    #     gridscale_spatial(
    #         d_gs_spatial,
    #         analysis,
    #         countries_mask,
    #         countries_regions,
    #         flags['extr'],
    #     )

    # checking fraction of countries emerged from noise (appears to decrease over GMT trajectories per birth year, which)
    # for step in GMT_labels:
    #     emerged_countries = xr.where(ds_ae_strj['age_emergence'].sel(GMT=step,birth_year=2000).mean(dim='run')>0,1,0).sum(dim='country') / len(ds_ae_strj.country.data)
    #     print(emerged_countries.item())
    
    

    # # rename ds_pop_frac['mean_unprec_all'] to 'unprec_Fraction

    # for row,var in zip(axes,[ds_le['lifetime_exposure'],ds_ae['age_emergence'],ds_pf['unprec_fraction']]):
    #     ax,cntry in zip(row,sample_countries):
    #         frame = {
                
    #         }
    #         df
            

    # collect all data arrays for age of emergence into dataset for finding age per birth year
    # ds_age_emergence = xr.merge([
    #     ds_age_emergence_15.rename({'age_emergence':'age_emergence_15'}),
    #     ds_age_emergence_20.rename({'age_emergence':'age_emergence_20'}),
    #     ds_age_emergence_NDC.rename({'age_emergence':'age_emergence_NDC'}),
    # ])
            
    # # plot pop frac of 3 main GMT mapped scenarios across birth years
    # plot_pop_frac_birth_year(
    #     ds_pop_frac_NDC,
    #     ds_pop_frac_15,
    #     ds_pop_frac_20,
    #     year_range,
    # )

    # # plot pop frac for 0.8-3.5 degree stylized trajectories across birth years
    # plot_pop_frac_birth_year_strj(
    #     ds_pop_frac_strj,
    #     df_GMT_strj,
    # )

    # plot pop frac and age emergence across GMT for stylized trajectories
    # top panel; (y: frac unprecedented, x: GMT anomaly @ 2100)
    # bottom panel; (y: age emergence, x: GMT anomaly @ 2100)
    # plots both approaches to frac unprecedented; exposed vs full cohorts
    # issue in these plots that early birth years for high warming trajectories won't live till 3-4 degree warming, but we misleadingly plot along these points
        # so, for e.g. 1970 BY, we need to limit the line up to life expectancy
        # will need cohort weighted mean of life expectancy across countries
        # 

    # # plot country mean age of emergence of 3 main GMT mapped scenarios across birth year
    # plot_age_emergence(
    #     ds_age_emergence_NDC,
    #     ds_age_emergence_15,
    #     ds_age_emergence_20,
    #     year_range,
    # )

    # # plot country mean age of emergence of stylized trajectories across birth years
    # plot_age_emergence_strj(
    #     ds_age_emergence_strj,
    #     df_GMT_strj,
    #     ds_cohorts,
    #     year_range,
    # )

    # calculate birth year emergence in simple approach
    # age emergences here shouldn't be from isolated 1.5, 2.0 and NDC runs; should have function that takes these main scenarios from stylized trajectories
    # gdf_exposure_emergence_birth_year = calc_exposure_emergence(
    #     ds_exposure,
    #     ds_exposure_pic,
    #     ds_age_emergence,
    #     gdf_country_borders,
    # )

    # # country-scale spatial plot of birth and year emergence
    # spatial_emergence_plot(
    #     gdf_exposure_emergence_birth_year,
    #     flags['extr'],
    #     flags['gmt'],
    # )

    # # plot stylized trajectories (GMT only)
    # plot_stylized_trajectories(
    #     df_GMT_strj,
    #     d_isimip_meta,
    #     year_range,
    # )

    # # plot pop frac across GMT for stylized trajectories; add points for 1.5, 2.0 and NDC from original analysis as test
    # plot_pop_frac_birth_year_GMT_strj_points(
    #     ds_pop_frac_strj,
    #     ds_age_emergence_strj,
    #     df_GMT_strj,
    #     ds_cohorts,
    #     ds_age_emergence,
    #     ds_pop_frac_15,
    #     ds_pop_frac_20,
    #     ds_pop_frac_NDC,
    #     ind_15,
    #     ind_20,
    #     ind_NDC,
    #     year_range,
    # )

# %%
