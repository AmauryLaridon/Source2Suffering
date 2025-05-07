# ---------------------------------------------------------------------------------------------------#
# Functions to compute and print the different reports needed                                        #
# ---------------------------------------------------------------------------------------------------#

#%%---------------------------------------------------------------#
# Libraries                                                       #
#-----------------------------------------------------------------#

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
scripts_dir, data_dir, data_dem4cli_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

#%% ---------------------------------------------------------------------#
# Sample analytics for Assessment                                        #
# -----------------------------------------------------------------------#

#%% ----------------------------------------------------------------#
# Generating large latex tables on CF data per country, birth year  #
# and 1.5, 2.5 and 3.5 degree scenario                              #
# ------------------------------------------------------------------#

def print_latex_table_unprecedented(
    flags,
    da_gs_popdenom,
):

    # input
    unprec_level="unprec_99.99"      
    bys=np.arange(1960,2021,10)
    # bys=np.arange(1960,2021,1)
    gmt_legend={
        GMT_indices_plot[0]:'1.5',
        GMT_indices_plot[1]:'2.5',
        GMT_indices_plot[2]:'3.5',
    }
    extremes = [
        'burntarea', 
        'cropfailedarea', 
        'driedarea', 
        'floodedarea', 
        'heatwavedarea', 
        'tropicalcyclonedarea',
    ]
    extremes_labels = {
        'burntarea': 'wildfires',
        'cropfailedarea': 'crop failures',
        'driedarea': 'droughts',
        'floodedarea': 'floods',
        'heatwavedarea': 'heatwaves',
        'tropicalcyclonedarea': 'tropical cyclones',
    }  

    for extr in extremes:
        
        # open dictionary of metadata for sim means and CF data per extreme
        with open(data_dir+'{}/{}/isimip_metadata_{}_{}_{}.pkl'.format(flags['version'],extr,extr,flags['gmt'],flags['rm']), 'rb') as file:
            d_isimip_meta = pk.load(file)     
        with open(data_dir+'{}/{}/gridscale_aggregated_pop_frac_{}.pkl'.format(flags['version'],extr,extr), 'rb') as f:
            ds_pf_gs = pk.load(f)           

        sims_per_step = {}
        sims_per_step[extr] = {}
        for step in GMT_indices_plot:
            sims_per_step[extr][step] = []
            for i in list(d_isimip_meta.keys()):
                if d_isimip_meta[i]['GMT_strj_valid'][step]:
                    sims_per_step[extr][step].append(i)      
        
        da_p_gs_plot = ds_pf_gs[unprec_level].loc[{
            'GMT':GMT_indices_plot,
        }]
        df_list_gs = []
        for step in GMT_indices_plot:
            da_p_gs_plot_step = da_p_gs_plot.loc[{'run':sims_per_step[extr][step],'GMT':step}].mean(dim='run')
            da_cf_gs_plot_step = da_p_gs_plot_step / da_gs_popdenom * 100
            df_cf_gs_plot_step = da_cf_gs_plot_step.to_dataframe(name='CF').reset_index()
            df_cf_gs_plot_step['P'] = da_p_gs_plot_step.to_dataframe(name='P').reset_index().loc[:,['P']] / 1000
            df_list_gs.append(df_cf_gs_plot_step)
        df_cf_gs_plot = pd.concat(df_list_gs)

        df_cf_gs_plot['CF'] = df_cf_gs_plot['CF'].fillna(0).round(decimals=0).astype('int') 
        df_cf_gs_plot['P'] = df_cf_gs_plot['P'].fillna(0).round(decimals=0).astype('int') 
        df_cf_gs_plot['P (CF)'] = df_cf_gs_plot.apply(lambda x: '{} ({})'.format(str(x.P),str(x.CF)), axis=1)

        # print latex per step
        for step in GMT_indices_plot:
            
            print('')
            print('Running latex print of CF for {} under {} pathway'.format(extremes_labels[extr],gmt_legend[step]))
            print('')
            
            df_latex = df_cf_gs_plot[df_cf_gs_plot['GMT']==step].copy()

            if flags['gridscale_country_subset']:

                df_latex = df_latex[df_latex["country"] == countries]
                df_cntry_by = df_latex.loc[:, ['country', 'birth_year', 'P (CF)']].set_index('country')
            else:
                df_cntry_by = df_latex.loc[:,['country','birth_year','P (CF)']].set_index('country')

            for by in bys:
                df_cntry_by[by] = df_cntry_by[df_cntry_by['birth_year']==by].loc[:,['P (CF)']]
                
            df_cntry_by = df_cntry_by.drop(columns=['birth_year','P (CF)']).drop_duplicates() 

            # latex
            caption = '\\caption{{\\textbf{{Absolute population (in thousands) of cohorts living unprecedented exposure to {0} and CF\\textsubscript{{{0}}} (\\%) per country and birth year in a {1}\\degree C pathway}}}}\\\\'.format(extremes_labels[extr],gmt_legend[step])
            headers = list(df_cntry_by.columns.astype('str'))
            headers = ['Country'] + headers
            data = {}
            for row in list(df_cntry_by.index):
                if len(str(row).split()) > 1:
                    newrow = ' \\\ '.join(str(row).split())
                    newrow = '\makecell[l]{{{}}}'.format(newrow)    
                    data[str(newrow)] = list(df_cntry_by.loc[row,:].values)
                else:
                    data[str(row)] = list(df_cntry_by.loc[row,:].values)

            textabular = f" l |{' c '*(len(headers)-1)}"
            texheader = " & ".join(headers) + "\\\\"
            texdata = "\\hline\n"

            for label in data:
                if label == "z":
                    texdata += "\\hline\n"
                texdata += f"{label} & {' & '.join(map(str,data[label]))} \\\\\n"

            print('\\small')
            print('\\begin{longtable}{'+textabular+'}')
            print(caption)
            print(texheader)
            print(texdata,end='')
            print('\\end{longtable}')
            print('\\normalsize') 
            print('\\clearpage')       