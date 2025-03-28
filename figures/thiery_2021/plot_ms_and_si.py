# -----------------------------------------------------------------------------------------------------------------
# Subscript to defines the functions needed for the plots from Thiery et al.(2021)
#
# This subscript is based on ms_plotscript.m from Thiery et al.(2021) is a translation of the main figures 
# in order to validate the translation of Thiery's Matlab framework in Python to obtain backward compatibility
# -----------------------------------------------------------------------------------------------------------------

#%%  ----------------------------------------------------------------
# Libraries 
# -------------------------------------------------------------------

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


print("Framework to reproduce the figures from Thiery et al. (2021) not yet developped in python.")

#%%  ----------------------------------------------------------------
# Initialization 
# -------------------------------------------------------------------

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

#%%  ----------------------------------------------------------------
# Plot f1
# -------------------------------------------------------------------

def plot_fig1():

    # Assuming you have already loaded your data in the respective variables (landfrac_peryear_15, landfrac_peryear_20, etc.)

    # Example data for demonstration (Replace with actual data)
    years_SR15 = np.arange(1960, 2100)
    ages = np.arange(0, 101)

    # Generate random values for land fraction per year
    landfrac_peryear_15 = np.random.random((5, 12, len(years_SR15))) * 0.5
    landfrac_peryear_20 = np.random.random((5, 12, len(years_SR15))) * 0.5
    landfrac_peryear_NDC = np.random.random((5, 12, len(years_SR15))) * 0.5

    # Define colors (similar to MATLAB setup)
    color_15 = 'blue'
    color_20 = 'green'
    color_NDC = 'red'
    boxcolor = 'black'

    # Smooth function for the plots (Gaussian filter)
    def smooth(data, window_size):
        return gaussian_filter1d(data, sigma=window_size)

    # Select region and extreme to analyze
    ind_extreme = 4  # Ensure this is within the range 0 to 4 for axis 0
    ind_region = 11  # Ensure this is within the range 0 to 11 for axis 1
    age_ref_plot = 60
    age_young = 20

    # Prepare data (panel a - smoothed time series)
    landfrac_15_plot = smooth(landfrac_peryear_15[ind_extreme, ind_region, :], 10) * 100
    landfrac_20_plot = smooth(landfrac_peryear_20[ind_extreme, ind_region, :], 10) * 100
    landfrac_NDC_plot = smooth(landfrac_peryear_NDC[ind_extreme, ind_region, :], 10) * 100

    # Panel a - Plot the time series with uncertainty bands
    fig, ax1 = plt.subplots(figsize=(10, 5))

    ax1.plot(years_SR15[:-4], landfrac_15_plot[:-4], color=color_15, linewidth=2)
    ax1.plot(years_SR15[:-4], landfrac_20_plot[:-4], color=color_20, linewidth=2)
    ax1.plot(years_SR15[:-4], landfrac_NDC_plot[:-4], color=color_NDC, linewidth=2)

    # Adding uncertainty bands (simulated with random data)
    landfrac_15_mms_plot = np.random.random(landfrac_15_plot.shape) * 0.1
    landfrac_20_mms_plot = np.random.random(landfrac_20_plot.shape) * 0.1
    landfrac_NDC_mms_plot = np.random.random(landfrac_NDC_plot.shape) * 0.1

    ax1.fill_between(years_SR15[:-4], landfrac_15_plot[:-4] - landfrac_15_mms_plot[:-4] / 2,
                    landfrac_15_plot[:-4] + landfrac_15_mms_plot[:-4] / 2, color=color_15, alpha=0.1)
    ax1.fill_between(years_SR15[:-4], landfrac_20_plot[:-4] - landfrac_20_mms_plot[:-4] / 2,
                    landfrac_20_plot[:-4] + landfrac_20_mms_plot[:-4] / 2, color=color_20, alpha=0.1)
    ax1.fill_between(years_SR15[:-4], landfrac_NDC_plot[:-4] - landfrac_NDC_mms_plot[:-4] / 2,
                    landfrac_NDC_plot[:-4] + landfrac_NDC_mms_plot[:-4] / 2, color=color_NDC, alpha=0.1)

    ax1.set_xlabel('Time')
    ax1.set_ylabel('Land area annually exposed to heatwaves [%]')
    ax1.set_xticks(np.arange(1960, 2100, 20))
    ax1.set_ylim(0, 55)

    # Panel b - Bar plot with exposure data
    exposure_bars = np.random.random((2, 3)) * 4  # Example bar data

    fig, ax2 = plt.subplots(figsize=(5, 5))
    categories = ['Age 60', 'Age 20']
    bar_width = 0.3
    x = np.arange(len(categories))  # Positions for bars

    # Plot bars for each scenario
    for i in range(exposure_bars.shape[1]):  # Loop over the number of scenarios (3)
        ax2.bar(x + i * bar_width, exposure_bars[:, i], width=bar_width, label=f'Scenario {i+1}')

    ax2.set_xticks(x + bar_width)
    ax2.set_xticklabels(categories)
    ax2.set_xlabel('Age of person in 2020')
    ax2.set_ylabel('Lifetime heatwave exposure [-]')
    ax2.legend()

    # Add EMF data (simulated for demonstration)
    for i in range(exposure_bars.shape[1]):
        for j in range(len(categories)):
            ax2.text(x[j] + i * bar_width, exposure_bars[j, i], f'x{np.random.randint(1, 20)}',
                    ha='center', va='bottom', color='black', fontsize=10)

    # Panel c - Plot EMF with boxplot
    fig, ax3 = plt.subplots(figsize=(5, 5))
    EMF_plot_15 = np.random.random(101) * 5  # Example EMF data
    EMF_plot_20 = np.random.random(101) * 5
    EMF_plot_NDC = np.random.random(101) * 5

    ax3.plot(ages, EMF_plot_15, color=color_15, linewidth=3)
    ax3.plot(ages, EMF_plot_20, color=color_20, linewidth=3)
    ax3.plot(ages, EMF_plot_NDC, color=color_NDC, linewidth=3)

    ax3.set_xlabel('Age of person in 2020')
    ax3.set_ylabel('Exposure multiplication factor over lifetime [-]')

    # Boxplot-like visualization
    ax3.boxplot([EMF_plot_15, EMF_plot_20, EMF_plot_NDC], positions=[1, 5, 9])

    plt.tight_layout()
    plt.savefig(scripts_dir+'/figures/thiery_2021/f1.png',dpi=1000,bbox_inches='tight')
    #plt.show()


