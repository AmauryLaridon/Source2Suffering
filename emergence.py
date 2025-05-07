# ----------------------------------------------------------------------------#
# Subscript to execute the emergence analysis                                 #
# ----------------------------------------------------------------------------#
      
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
import mapclassify as mc
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import interpolate
import cartopy.crs as ccrs
from settings import *
scripts_dir, data_dir, data_dem4clim_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()

global grid_area
grid_area = xr.open_dataarray(data_dir+'isimip/grid_resolution/clm45_area.nc4')

sys.path.append(os.path.abspath(scripts_dir+"/pf_scripts"))
from pf_emergence import *

# ----------------------------------------------------------------#
# Execute emergence                                               #
# ----------------------------------------------------------------#

if flags['emergence']:

    # --------------------------------------------------------------------
    # process emergence of cumulative exposures, mask cohort exposures for time steps of emergence

    print("Computing Emergence of cumulative exposures")
    
    if flags['birthyear_emergence']:
        
        by_emergence = np.arange(1960,2101)
        
    else:
        
        by_emergence = birth_years        
    
    if not os.path.isfile(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version'])):
        
        # need new cohort dataset that has total population per birth year (using life expectancy info; each country has a different end point)
        da_cohort_aligned = calc_birthyear_align(
            da_cohort_size,
            df_life_expectancy_5,
            by_emergence,
        )
        
        # convert to dataset and add weights
        ds_cohorts = ds_cohort_align(
            da_cohort_size,
            da_cohort_aligned,
        )
        
        # pickle birth year aligned cohort sizes and global mean life expectancy
        with open(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version']), 'wb') as f:
            pk.dump(ds_cohorts,f)  

    else:
        
        # load pickled birth year aligned cohort sizes and global mean life expectancy
        with open(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version']), 'rb') as f:
            ds_cohorts = pk.load(f)                             
    
    ds_ae_strj, ds_pf_strj = strj_emergence(
        d_isimip_meta,
        df_life_expectancy_5,
        ds_exposure_pic,
        ds_cohorts,
        by_emergence,
        flags,
    )
        
else: # load pickles

    print("Loading Emergence of cumulative exposures")
    
    pass
    
    # birth year aligned population
    with open(data_dir+'{}/country/cohort_sizes.pkl'.format(flags['version']), 'rb') as f:
        ds_cohorts = pk.load(f)
    
    # pop frac
    with open(data_dir+'{}/{}/pop_frac_{}.pkl'.format(flags['version'],flags['extr'],flags['extr']), 'rb') as f:
        ds_pf_strj = pk.load(f)         