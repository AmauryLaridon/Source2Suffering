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
import regionmask as rm
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import interpolate
import cartopy.crs as ccrs
import shutil
import glob
from settings import *

ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins = init()

#%% ----------------------------------------------------------------
# grid scale
# ------------------------------------------------------------------

def add_directories(
    list_countries,
    flags,
):
    
    extremes = [
        'burntarea', 
        'cropfailedarea', 
        'driedarea', 
        'floodedarea', 
        'heatwavedarea', 
        'tropicalcyclonedarea',
    ]
    
    for extr in extremes:

        for cntry in list_countries:
            
            if not os.path.exists('./data/{}/{}/{}'.format(flags['version'],extr,cntry)):
                os.makedirs('./data/{}/{}/{}'.format(flags['version'],extr,cntry)) # testing makedirs
                
            for file in glob.glob('./data/{}/{}/gridscale_le_pic_{}_{}.pkl'.format(flags['version'],extr,extr,cntry)):
                shutil.move(file, './data/{}/{}/{}'.format(flags['version'],extr,cntry))            
                    
            for file in glob.glob('./data/{}/{}/gridscale_pic_qntls_{}_{}.pkl'.format(flags['version'],extr,extr,cntry)):      
                shutil.move(file, './data/{}/{}/{}'.format(flags['version'],extr,cntry))          
                                                
            # check for pickle of gridscale lifetime exposure (da_le); process if not existing; os.mkdir('./data/{}/{}/{}'.format(flags['version'],extr,cntry))
            if extr != 'heatwavedarea':
                for file in glob.glob('./data/{}/{}/gridscale_le_{}_{}_*.pkl'.format(flags['version'],extr,extr,cntry)):
                    shutil.move(file, './data/{}/{}/{}'.format(flags['version'],extr,cntry))
                
            for file in glob.glob('./data/{}/{}/gridscale_le_{}_{}_*.pkl'.format(flags['version'],extr,extr,cntry)):
                shutil.move(file, './data/{}/{}/{}'.format(flags['version'],extr,cntry))                 
                            
            for file in glob.glob('./data/{}/{}/gridscale_emergence_mask_{}_{}_*.pkl'.format(flags['version'],extr,extr,cntry)):
                shutil.move(file, './data/{}/{}/{}'.format(flags['version'],extr,cntry))  
