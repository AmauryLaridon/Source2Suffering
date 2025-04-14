# -----------------------------------------------------------------------------------------------------------------#
# Functions needed for plots in the Source2Suffering project                                                       #
# -----------------------------------------------------------------------------------------------------------------#

#%%  ----------------------------------------------------------------#
# Libraries                                                          #
# -------------------------------------------------------------------#

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
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import regionmask as rm
import geopandas as gpd
from scipy import interpolate
from scipy import stats as sts
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import uniform_filter1d
import cartopy.crs as ccrs
import seaborn as sns
import cartopy as cr
import cartopy.feature as feature
import imageio
#from skimage import exposure

from settings import *
scripts_dir, data_dir, ages, age_young, age_ref, age_range, year_ref, year_start, birth_years, year_end, year_range, GMT_max, GMT_min, GMT_inc, RCP2GMT_maxdiff_threshold, year_start_GMT_ref, year_end_GMT_ref, scen_thresholds, GMT_labels, GMT_window, GMT_current_policies, pic_life_extent, nboots, resample_dim, pic_by, pic_qntl, pic_qntl_list, pic_qntl_labels, sample_birth_years, sample_countries, GMT_indices_plot, birth_years_plot, letters, basins, countries = init()


#%%------------------------------------------------------------------#
# Initialization                                                     #
# -------------------------------------------------------------------#

# Set color scale axes
caxes = {
    "EMF": [1/5, 5],  # Used for geometric mean - revised submission
    "EMF_heatwaves": [1/10, 10],  # Used for heatwaves - revised submission
    "EMF_BE": [1, 4],
    "PCT_BE_plot": [0.001, 100000],
    "EMF_BE_coldwaves": [1/10, 1],
    "dexposure": [-20, 20],
    "regions_geographic": [0.5, 7.5],
    "regions_income": [0.5, 4.5]
}

# Set colormaps
colormaps = {
    "EMF": plt.get_cmap("coolwarm", 28),  # Custom colormap selection
    "EMF_BE_coldwaves": plt.get_cmap("Blues", 40),
    "regions_geographic": plt.get_cmap("tab10"),
    "regions_income": plt.get_cmap("Pastel1"),
    "dexposure": plt.get_cmap("coolwarm", 20),
    "pie": plt.get_cmap("Set2")
}

# Adjust colormap to avoid close resemblance of dark red
colormap_array = colormaps["EMF"](np.linspace(0, 1, 28))
colormap_array[0] = colormap_array[1]
colormaps["EMF"] = mcolors.ListedColormap(colormap_array)

# Define custom BE colormap
colormaps["BE_4"] = np.array([
    [255/255, 255/255, 255/255],  # White
    [252/255, 204/255, 18/255],  # Yellow
    [218/255, 60/255, 48/255],  # Red
    [143/255, 34/255, 91/255]  # Dark Purple
])

# Fonction to rescale the intensities (similar to exposure.rescale_intensity)
def rescale_intensity(arr, out_range=(0, 1)):
    min_val = np.min(arr)
    max_val = np.max(arr)
    return (arr - min_val) / (max_val - min_val) * (out_range[1] - out_range[0]) + out_range[0]


colormaps["BE_pseudodiscrete"] = np.repeat(colormaps["BE_4"], 20, axis=0)
colormaps["BE_pseudodiscrete"] = rescale_intensity(colormaps["BE_pseudodiscrete"], out_range=(0, 1))

colormaps["BE_continuous"] = rescale_intensity(colormaps["BE_4"], out_range=(0, 1))

# Define color settings
colors = plt.get_cmap("tab20").colors
color_15, color_20, color_NDC = colors[5], colors[6], colors[7]

# Define axes and sea color
darkcolor = [0.3, 0.3, 0.3]  # Dark grey
axcolor = [0.5, 0.5, 0.5]  # Normal grey
boxcolor = [0.7, 0.7, 0.7]  # Light grey
seacolor = [0.95, 0.95, 0.95]  # Very light grey

# Define alphabet letters
alphabet = [chr(i) for i in range(ord('a'), ord('z')+1)]

# Define panel letters for region plots
panelletters_regions = ['b', 'c', 'd', 'd', 'a', 'b', 'e', 'f', 'g', 'h', 'c', 'a']

# Define pictogram dimensions
pictogram_dims = np.array([
    [-2.8, -2.2, -0.2, -1.0],  # Burnt area
    [-3.0, -2.0, -0.2, -1.2],  # Crop failed area
    [-2.8, -2.2, -0.2, -1.0],  # Dried area
    [-2.8, -2.0, -0.2, -1.0],  # Flooded area
    [-2.8, -2.0, -0.2, -1.0],  # Heatwaved area
    [-2.8, -2.0, -0.2, -1.0],  # Tropical cyclone area
    [-3.0, -1.7, -0.2, -1.0]   # All events
])

# Function to adjust gamma
def adjust_gamma(image, gamma=1.0):
    return np.power(image, gamma)

# Load and adjust pictograms
pictogram_all_2x3 = imageio.imread(scripts_dir+'/figures/thiery_2021/pictograms/pictogram_all_2x3.png')
pictogram_all_2x3 = adjust_gamma(1 - pictogram_all_2x3 / 255, gamma=0.7)

pictogram_all_1x6 = imageio.imread(scripts_dir+'/figures/thiery_2021/pictograms/pictogram_all_1x6.png')
pictogram_all_1x6 = adjust_gamma(1 - pictogram_all_1x6 / 255, gamma=0.7)

# Define world regions and income group indices
ind_worldreg = [8, 4, 2, 7, 10, 1, 9, 12]  # World regions from West to East
ind_income = [5, 6, 11, 3]  # Income categories from low to high

# Define line colors
fig, ax = plt.subplots()
colormaps["worldregions"] = plt.get_cmap("Dark2", 8)
linecolors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
linecolors.append(colormaps["worldregions"](7))
plt.close(fig)


#%%------------------------------------------------------------------#
# Functions                                                          #
# -------------------------------------------------------------------#

def plot_dev_fig1(flags, ds, country, run, GMT):
    """
    Plot lifetime exposure for a specific country, run, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - run (int): Numéro de simulation (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """
    
    exposure = ds['lifetime_exposure'].sel(
        country=country,
        run=run,
        GMT=GMT  
    )

    
    birth_years = ds['birth_year'].values

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.title("{} : Lifetime Exposure in {} (Run {}, GMT {})".format(flags['extr'], country, run, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/fig1_lifetime_exposure_{}_GMT_{}_isimip_sim{}.png'.format(country,GMT,run))

def plot_dev_fig2(ds, var_name, coords_dict):
    """
    Function to plot a variable from the dataset ds based on given coordinates.

    Parameters:
    - ds (xr.Dataset): The dataset containing the variable to plot.
    - var_name (str): The name of the variable to plot (e.g., 'exposure_trend_ar6', 'mean_exposure_trend_ar6', etc.).
    - coords_dict (dict): Dictionary of coordinates for slicing the dataset. The dictionary should contain 
                          keys corresponding to the coordinate names and their respective values.
                          Example: {'run': 1, 'GMT': 0, 'region': 'Africa', 'year': 2000}

    Returns:
    - None
    """

    # Ensure the variable exists in the dataset
    if var_name not in ds.data_vars:
        print(f"Error: The variable '{var_name}' does not exist in the dataset.")
        return

    # Apply the slicing based on the provided coordinates
    sliced_data = ds[var_name]

    # Create the title showing the coordinates being used for slicing
    title = f"Plot of {var_name} for " + ", ".join([f"{k}={v}" for k, v in coords_dict.items()])
    
    # Loop over the coordinates in the dictionary and apply them
    for coord, value in coords_dict.items():
        if coord not in sliced_data.coords:
            print(f"Error: The coordinate '{coord}' does not exist in the variable's dimensions.")
            return
        # Slice the data for each coordinate
        sliced_data = sliced_data.sel({coord: value}, drop=True)

    # Plot the data (Assuming the remaining dimensions after slicing are 'year' and 'value')
    plt.figure(figsize=(10, 6))
    plt.plot(sliced_data['year'], sliced_data.values)
    plt.xlabel('Year')
    plt.ylabel(var_name)
    plt.title(title)
    plt.grid(True)
    # Save the plot with a dynamic filename based on the coordinates
    if var_name == 'exposure_trend_ar6':
        filename = scripts_dir + '/figures/source2suffering/development/fig2_lifetime_exposure_trends_isimip_sim{}_GMT_{}_region_{}.png'.format(
        coords_dict.get('run', 'NA'),coords_dict.get('GMT', 'NA'), coords_dict.get('region', 'NA'))
    plt.savefig(filename)

def plot_dev_fig3(ds_regions, flags, ds, region, run, GMT):
    """
    Plot lifetime exposure for a specific region, run, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - run (int): Numéro de simulation (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """
    
    exposure = ds['lifetime_exposure'].sel(
        region=region,
        run=run,
        GMT=GMT  
    )

    
    birth_years = ds['birth_year'].values
    region_name = ds_regions['name'].sel(region=region)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.title("{} : Lifetime Exposure in {} (Run {}, GMT {})".format(flags['extr'], region_name, run, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/fig3_lifetime_exposure_region_{}_GMT_{}_isimip_sim{}.png'.format(region,GMT,run))