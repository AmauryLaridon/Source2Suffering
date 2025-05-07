# ---------------------------------------------------------------------------------------------#
# Functions to compute and print the different reports needed in the Source2Suffering Project  #
# ---------------------------------------------------------------------------------------------#

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

#%% ----------------------------------------------------------------#
# Sample analytics for paper                                        #
# ------------------------------------------------------------------#

print("Framework to reproduce the reports for the Source2Suffering Project not yet developped in python.")
