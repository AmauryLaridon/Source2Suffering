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
import imageio.v2 as imageio

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

def plot_dev_fig1(flags, extr, ds, country, run, GMT):
    """
    Plot lifetime exposure for a specific country, run, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - run (int): Numéro de simulation (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """

    plt.close('all') 
    
    exposure = ds['le_percountry_perrun_BE'].sel(
        country=country,
        run=run,
        GMT=GMT  
    )

    
    birth_years = ds['birth_year'].values

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.title("{} : Lifetime Exposure in {} (Run {}, GMT {})".format(extr, country, run, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig1_lifetime_exposure_{}_GMT_{}_isimip_sim_{}_gmt_{}_{}.png'.format(extr,country,GMT,run,flags['gmt'],flags['rm']))

def plot_dev_fig2(flags, ds, var_name, coords_dict):
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

    plt.close('all') 

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
        filename = scripts_dir + '/figures/source2suffering/development/{}/fig2_lifetime_exposure_trends_isimip_sim_{}_GMT_{}_region_{}.png'.format(flags['extr'],
        coords_dict.get('run', 'NA'),coords_dict.get('GMT', 'NA'), coords_dict.get('region', 'NA'))
    plt.savefig(filename)


def plot_dev_fig3(ds_regions, flags, extr, ds, region, run, GMT):
    """
    Plot lifetime exposure for a specific region, run, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - run (int): Numéro de simulation (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """
    
    plt.close('all') 

    exposure = ds['le_perregion_perrun_BE'].sel(
        region=region,
        run=run,
        GMT=GMT  
    )

    
    birth_years = ds['birth_year'].values
    region_name = ds_regions['name'].sel(region=region)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.title("{} : Lifetime Exposure in {} (Run {}, GMT {})".format(extr, region_name, run, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig3_lifetime_exposure_region_{}_GMT_{}_isimip_sim_{}_gmt_{}_{}.png'.format(extr,region_name,GMT,run,flags['gmt'],flags['rm']))


def plot_dev_fig4(flags, extr, ds, country, GMT):
    """
    Plot MMM lifetime exposure for a specific country, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - GMT (str or int): Niveau GMT (ex: '21')
    """

    plt.close('all') 
    
    exposure = ds['mmm_BE'].sel(
        country=country,
        GMT=GMT  
    )

    y_std = ds['std_BE'].sel(
        country = country,
        GMT = GMT
    )/2

    birth_years = ds['birth_year'].values

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$')
    plt.title("{} : MMM Lifetime Exposure in {} (GMT {})".format(extr, country, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig4_mmm_lifetime_exposure_{}_GMT_{}_gmt_{}_{}.png'.format(extr,country,GMT,flags['gmt'],flags['rm']))


def plot_dev_fig5(ds_regions, flags, extr, ds, region, GMT):
    """
    Plot MMM lifetime exposure for a specific region, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - region (int): Region index (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """
    
    plt.close('all') 

    exposure = ds['mmm_BE'].sel(
        region = region,
        GMT = GMT  
    )

    birth_years = ds['birth_year'].values
    region_name = ds_regions['name'].sel(region=region)
    y_std = ds['std_BE'].sel(
        region = region,
        GMT = GMT
    )/2

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$')
    plt.title("{} : MMM Lifetime Exposure in {} (GMT {})".format(extr, region_name, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig5_mmm_lifetime_exposure_region_{}_GMT_{}_gmt_{}_{}.png'.format(extr,region_name,GMT,flags['gmt'],flags['rm']))

def plot_dev_fig6_country(flags, extr, ds, country, EMF):
    """
    Plot MMM lifetime exposure to heatwave for all countries, and three GMT values.

    Parameters:
    - ds (xarray.Dataset): Dataset containing the variable 'lifetime_exposure'
    - country (string): Name of the country (ex: 'Belgium')
    - EMF (Boolean): Choose to express the results in terms of EMF or remain in absolute lifetime exposure
    """

    plt.close('all') 

    plt.figure(figsize=(12, 8))

    if extr == 'burntarea':
        extr_name = 'Wildfires'
    if extr == 'cropfailedarea':
        extr_name = 'Crop failures'
    if extr == 'driedarea':
        extr_name = 'Droughts'
    if extr == 'floodedarea':
        extr_name = 'River floods'
    if extr == 'heatwavedarea':
        extr_name = 'Heatwaves'
    if extr == 'tropicalcyclonedarea':
        extr_name = 'Tropical Cyclones'

    if flags['gmt'] == 'ar6_new':
        GMT_list = [0, 11, 20]
        GMT_label = ['1.5°C', '2.5°C', '3.5°C']
    if flags['gmt'] == 'original':
        GMT_label = ['1.5°C', '2.0°C', 'NDC']

    GMT_color = ['tab:blue', 'tab:orange', 'tab:red']
    i = 0

    if EMF:

        for GMT in GMT_list:
            exposure = ds.sel(
                country=country,
                GMT=GMT
            )

            plt.plot(birth_years, exposure, linestyle='-', color=GMT_color[i], label=GMT_label[i])
            i += 1

        plt.title("{} : MMM Lifetime Exposure EMF in {}".format(extr, country))
        plt.xlabel("Birth Year")
        plt.ylabel("Exposure Multiplication Factor (EMF)")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(scripts_dir + '/figures/source2suffering/development/{}/storyline/fig6_mmm_EMF_{}_gmt_{}_{}.png'.format(extr, country, flags['gmt'], flags['rm']))

    else:

        if flags['gmt'] == 'ar6_new':

            exposure_data = []  # To store data for CSV export if needed

            for GMT in GMT_list:

                exposure = ds['mmm_BE'].sel(
                    country=country,
                    GMT=GMT
                ).values

                y_std = (ds['std_BE'].sel(
                    country=country,
                    GMT=GMT
                ).values) / 2

                plt.plot(birth_years, exposure, linestyle='-', color=GMT_color[i], label=GMT_label[i], lw=3)
                plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, color=GMT_color[i], lw=3)

                exposure_data.append((GMT_label[i], exposure, y_std))

                i += 1

        if flags['gmt'] == 'original':

            exposure_data = []

            exposure_15 = ds['mmm_15'].sel(country=country).values
            y_std_15 = ds['std_15'].sel(country=country).values / 2

            exposure_20 = ds['mmm_20'].sel(country=country).values
            y_std_20 = ds['std_20'].sel(country=country).values / 2

            exposure_NDC = ds['mmm_NDC'].sel(country=country).values
            y_std_NDC = ds['std_NDC'].sel(country=country).values / 2

            plt.plot(birth_years, exposure_15, linestyle='-', color=GMT_color[0], label=GMT_label[0], lw=3)
            plt.fill_between(birth_years, exposure_15 - y_std_15, exposure_15 + y_std_15, alpha=0.3, color=GMT_color[0], lw=3)

            plt.plot(birth_years, exposure_20, linestyle='-', color=GMT_color[1], label=GMT_label[1], lw=3)
            plt.fill_between(birth_years, exposure_20 - y_std_20, exposure_20 + y_std_20, alpha=0.3, color=GMT_color[1], lw=3)

            plt.plot(birth_years, exposure_NDC, linestyle='-', color=GMT_color[2], label=GMT_label[2], lw=3)
            plt.fill_between(birth_years, exposure_NDC - y_std_NDC, exposure_NDC + y_std_NDC, alpha=0.3, color=GMT_color[2], lw=3)

            exposure_data.append((GMT_label[0], exposure_15, y_std_15))
            exposure_data.append((GMT_label[1], exposure_20, y_std_20))
            exposure_data.append((GMT_label[2], exposure_NDC, y_std_NDC))

        plt.title("Lifetime Exposure to {}\nin {}".format(extr_name, country), fontsize=18, fontweight='bold')
        plt.xlabel("Birth Year", fontsize=16)
        plt.ylabel("Lifetime Exposure", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(True)
        plt.tight_layout()
        plt.legend(fontsize=16)
        plt.savefig(scripts_dir + '/figures/source2suffering/development/{}/storyline/fig6_mmm_lifetime_exposure_{}_gmt_{}_{}.png'.format(extr, country, flags['gmt'], flags['rm']))
        print("Plot performed for {}".format(country))

        # CSV export for Japan
        if (country == 'Japan') and (flags['gmt'] == 'ar6_new'):
            print('Saving .csv file with data of the figure for Nikkei Japan')

            csv_filename = scripts_dir + '/output/assessment/Nikkei Japan/Lifetime Exposure/{}_fig1.csv'.format(extr_name)

            with open(csv_filename, 'w') as f:
                # Write header
                header = 'Birth Year'
                for label in GMT_label:
                    header += f', Exposure {label}, Std.dev {label}'
                f.write(header + '\n')

                # Write data rows
                for idx, year in enumerate(birth_years):
                    row = f"{year}"
                    for (_, exposure_vals, std_vals) in exposure_data:
                        row += f", {exposure_vals[idx]}, {std_vals[idx]}"
                    f.write(row + '\n')

            print(f'Data exported to: {csv_filename}')



def plot_dev_fig6_region(ds_regions, flags, extr, ds, region, EMF):
    """
    Plot MMM lifetime exposure to heatwave for all regions, and three GMT values.

    Parameters:
    - ds_regions (xarray.Dataset): Dataset with all the countries per region
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - region (int): Index of the region (ex: '1')
    """
    
    plt.close('all') 

    plt.figure(figsize=(12, 8))

    if extr=='burntarea':
        extr_name = 'Wildfires'
    if extr=='cropfailedarea':
        extr_name = 'Crop failures'
    if extr=='driedarea':
        extr_name = 'Droughts'
    if extr=='floodedarea':
        extr_name = 'River floods'
    if extr=='heatwavedarea':
        extr_name = 'Heatwaves'
    if extr=='tropicalcyclonedarea':
        extr_name = 'Tropical Cyclones'

    if flags['gmt'] == 'ar6_new':

        GMT_list = [0,11,20]
        GMT_label = ['1.5°C', '2.5°C', '3.5°C']

    if flags['gmt'] == 'original':
    
        GMT_label = ['1.5°C', '2.0°C', 'NDC']
        
    GMT_color = ['tab:blue','tab:orange','tab:red']
    region_name = ds_regions['name'].sel(region=region)
    i = 0
    
    if EMF:

        for GMT in GMT_list:
    
            exposure = ds.sel(
                region=region,
                GMT=GMT  
            )

            # y_std = ds['std'].sel(
            #     region=region,
            #     GMT = GMT
            # )

            plt.plot(birth_years, exposure, linestyle='-',color=GMT_color[i],label=GMT_label[i])
            #plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i])

            i +=1

        plt.title("{} : MMM Lifetime Exposure EMF in {}".format(extr, region_name))
        plt.xlabel("Birth Year")
        plt.ylabel("Exposure Multiplication Factor (EMF)")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/storyline/fig6_mmm_EMF_region_{}_gmt_{}_{}.png'.format(extr,region_name,flags['gmt'],flags['rm']))

    else: 

        if flags['gmt']=='ar6_new':

            for GMT in GMT_list:

                exposure = ds['mmm_BE'].sel(
                    region=region,
                    GMT=GMT  
                )

                y_std = ds['std_BE'].sel(
                    region=region,
                    GMT = GMT
                )/2

                plt.plot(birth_years, exposure, linestyle='-',color=GMT_color[i],label=GMT_label[i],lw=3)
                plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i],lw=3)

                i +=1
        
        if flags['gmt']=='original':

            exposure_15 = ds['mmm_15'].sel(
                    region=region,
                )

            y_std_15 = ds['std_15'].sel(
                region=region,
            )/2

            exposure_20 = ds['mmm_20'].sel(
                region=region,
                )

            y_std_20 = ds['std_20'].sel(
                region=region,
            )/2

            exposure_NDC = ds['mmm_NDC'].sel(
                region=region,
                )

            y_std_NDC = ds['std_NDC'].sel(
                region=region,
            )/2

            plt.plot(birth_years, exposure_15, linestyle='-',color=GMT_color[0],label=GMT_label[0], lw=3)
            plt.fill_between(birth_years, exposure_15 - y_std_15, exposure_15 + y_std_15, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[0], lw=3)

            plt.plot(birth_years, exposure_20, linestyle='-',color=GMT_color[1],label=GMT_label[1], lw=3)
            plt.fill_between(birth_years, exposure_20 - y_std_20, exposure_20 + y_std_20, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[1], lw=3)

            plt.plot(birth_years, exposure_NDC, linestyle='-',color=GMT_color[2],label=GMT_label[2], lw=3)
            plt.fill_between(birth_years, exposure_NDC - y_std_NDC, exposure_NDC + y_std_NDC, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[2], lw=3)

        plt.title("Lifetime Exposure to {} \n in the {} region".format(extr_name, region_name),fontsize=18,fontweight='bold')
        plt.xlabel("Birth Year",fontsize=16)
        plt.ylabel("Lifetime Exposure",fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(True)
        plt.tight_layout()
        plt.legend(fontsize=16)
        plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/storyline/fig6_mmm_lifetime_exposure_region_{}_gmt_{}_{}.png'.format(extr,region_name,flags['gmt'],flags['rm']))
        print("Plot performed for the {} region".format(region_name))

def plot_dev_fig7(
    ds_regions,
    df_GMT_strj,
    ds_le_perregion,
    extr,
    flags,
    region_ind
):
    """
    Plot one Burning Ember (BE) of lifetime exposure for heatwaves per birth years for a specific region.
    The variable as output for the z-axis is the lifetime exposure value. 

    Parameters:
    - ds_regions (DataSet) : contains all the informations about the regions 
    - df_GMT_strj (DataFrame) : Contains the GMT stylized trajectories 
    - flags (Array) : contains the parameters of the configuration
    - region_ind (int) : Index of the region of interest
    """

    plt.close('all') 


    region_name = ds_regions['name'].sel(region=region_ind)

    # labels for GMT ticks

    if flags['gmt']=='ar6_new':

        GMT_indices_ticks=[0,5,10,15,20]

        with open(data_dir+'temperature_trajectories_AR6/df_GMT_strj.pkl', 'rb') as f:
            df_GMT_plots = pk.load(f)


    if flags['gmt']=='original':

        GMT_indices_ticks=[7,12,17,24,27]

        with open(data_dir+'temperature_trajectories_SR15/df_GMT_strj.pkl', 'rb') as f:
            df_GMT_plots = pk.load(f)

    gmts2100 = np.round(df_GMT_plots.loc[2100,GMT_indices_ticks].values,1)    

 
    #------------------------------ Plot ----------------------------------#

    # Ticks style

    mpl.rcParams['xtick.labelcolor'] = 'black'
    mpl.rcParams['ytick.labelcolor'] = 'black'

    # Size of the figure

    fig = plt.figure(figsize=(8, 7))

    # Label of the hazard
    hazard_labels = {
        'burntarea': 'Wildfires',
        'cropfailedarea': 'Crop failures',
        'driedarea': 'Droughts',
        'floodedarea': 'River floods',
        'heatwavedarea': 'Heatwaves',
        'tropicalcyclonedarea': 'Tropical Cyclones'
    }
    extr_name = hazard_labels.get(extr, extr)

    # Init one axe 

    ax = plt.gca()

    # Main plot

    da = ds_le_perregion['mmm_BE'].loc[{
        'region':region_ind,
        'birth_year': np.arange(1960, 2021),
    }]

    max_val = da.max(skipna=True).item()

    if np.isnan(max_val):

        print("Plot not performed because NaN values are present in Lifetime Exposure for {}".format(region_name))

    else: 

        if extr=='heatwavedarea':

            max_val = da.max(skipna=True).item()

            levels_hw=np.arange(0,max_val,0.1)

            p =da.plot.contourf(
                x='birth_year',
                y='GMT',
                ax=ax,
                add_labels=False,
                levels=levels_hw,
                cmap='Reds',
                cbar_kwargs={'ticks': np.arange(0, max_val, 2)}
            )

        else:

            # Compute the min and max values of the data
            min_val = da.min(skipna=True).item()
            max_val = da.max(skipna=True).item()
            range_val = max_val - min_val

            # Early exit if data has no range
            if range_val == 0:
                print("Plot not performed: constant data for {}".format(region_name))
                return

            # --- Contour levels ---
            target_n_levels = 100
            max_n_levels = 200

            step = range_val / target_n_levels if target_n_levels != 0 else 0.01
            min_step = 0.001
            rounded_step = max(np.round(step, 3), min_step)

            n_levels = int(np.ceil(range_val / rounded_step)) + 1

            if n_levels > max_n_levels:
                rounded_step = range_val / max_n_levels
                n_levels = max_n_levels

            levels_hw = np.arange(min_val, min_val + rounded_step * n_levels, rounded_step)

            # --- Colorbar ticks (aligned with contour levels) ---
            target_n_ticks = 8
            tick_indices = np.linspace(0, len(levels_hw) - 1, target_n_ticks).astype(int)
            ticks_z = np.round(levels_hw[tick_indices], 2)

            cbar_kwargs = {'ticks': ticks_z}

            # --- Plot ---
            p = da.plot.contourf(
                x='birth_year',
                y='GMT',
                ax=ax,
                add_labels=False,
                levels=levels_hw,
                cmap='Reds',
                extend='both',
                cbar_kwargs=cbar_kwargs
            )

        # Set ticks and labels
        ax.set_yticks(
            ticks=GMT_indices_ticks,
            labels=gmts2100,
            color='black',
        )

        ax.set_xticks(
            ticks=np.arange(1960, 2025, 10),
            color='black',
        )  

        ax.tick_params(axis='both', labelsize=14)  

        # Title and aesthetics
        ax.set_title(
        "Lifetime Exposure to {} \n in the {} region".format(extr_name, region_name),
        fontsize=18, fontweight='bold', loc='center', color='black'
        )

        ax.spines['right'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')  

        ax.set_xlabel('Birth year', fontsize=16, color='black')
        ax.set_ylabel('GMT change relative to PI in 2100 [°C]', fontsize=16, color='black', rotation='vertical')

        plt.tight_layout()

        if flags['old_demo']:

            plt.savefig(
                scripts_dir + '/figures/source2suffering/sandbox/old_demography/fig7_BE_{}_lifetime_exposure_region_{}_gmt_{}_{}.png'.format(extr, region_name, flags['gmt'],flags['rm']
                )
            )

        else: 

            plt.savefig(
                scripts_dir + '/figures/source2suffering/development/{}/fig7_BE_{}_lifetime_exposure_region_{}_gmt_{}_{}.png'.format(
                    extr, extr, region_name, flags['gmt'],flags['rm']
                )
            )

    print("Plot performed for {}".format(region_name))


def plot_dev_fig8(
    df_GMT_strj,
    ds_le_percountry,
    flags,
    extr, 
    country
):
    """
    Plot one Burning Ember (BE) of lifetime exposure for heatwaves per birth years for a specific country.
    The variable as output for the z-axis is the lifetime exposure value. 

    Parameters:
    - df_GMT_strj (DataFrame) : Contains the GMT stylized trajectories 
    - flags (Array) : contains the parameters of the configuration
    - country (String) : Name of the country of interest
    """

    plt.close('all') 


    # labels for GMT ticks

    if flags['gmt']=='ar6_new':

        GMT_indices_ticks=[0,5,10,15,20]

        with open(data_dir+'temperature_trajectories_AR6/df_GMT_strj.pkl'.format(flags['rm']), 'rb') as f:
            df_GMT_plots = pk.load(f)


    if flags['gmt']=='original':

        GMT_indices_ticks=[7,12,17,24,27]

        with open(data_dir+'temperature_trajectories_SR15/df_GMT_strj.pkl'.format(flags['rm']), 'rb') as f:
            df_GMT_plots = pk.load(f)

    gmts2100 = np.round(df_GMT_plots.loc[2100,GMT_indices_ticks].values,1)     
 
    #------------------------------ Plot ----------------------------------#

    # Ticks style

    mpl.rcParams['xtick.labelcolor'] = 'black'
    mpl.rcParams['ytick.labelcolor'] = 'black'

    # Size of the figure

    fig = plt.figure(figsize=(8, 7))

    # Label of the hazard
    hazard_labels = {
        'burntarea': 'Wildfires',
        'cropfailedarea': 'Crop failures',
        'driedarea': 'Droughts',
        'floodedarea': 'River floods',
        'heatwavedarea': 'Heatwaves',
        'tropicalcyclonedarea': 'Tropical Cyclones'
    }
    extr_name = hazard_labels.get(extr, extr)

    # Init one axe 

    ax = plt.gca()

    # Main plot

    da = ds_le_percountry['mmm_BE'].loc[{
        'country':country,
        'birth_year': np.arange(1960, 2021),
    }]

    max_val = da.max(skipna=True).item()

    if np.isnan(max_val):

        print("Plot not performed because NaN values are present in Lifetime Exposure for {}".format(country))
    
    else:

        if extr=='heatwavedarea':
        
            levels_hw=np.arange(0,max_val,0.1)

            p =da.plot.contourf(
                x='birth_year',
                y='GMT',
                ax=ax,
                add_labels=False,
                levels=levels_hw,
                cmap='Reds',
                cbar_kwargs={'ticks': np.arange(0, max_val, 2)}
            )
        
        else:
    
            # Compute the min and max values of the data
            min_val = float(da.min(skipna=True))
            max_val = float(da.max(skipna=True))
            range_val = max_val - min_val

            # Early exit if data has no range
            if np.isclose(range_val, 0):
                print(f"Plot not performed: constant data for {country}")
                return

            # Define desired number of levels and ticks
            target_n_levels = 100
            max_n_levels = 200
            target_n_ticks = 10
            max_n_ticks = 20

            # Compute step size for contour levels
            step = range_val / target_n_levels
            min_step = 0.001  # enforce some visual distinction
            rounded_step = max(np.round(step, 3), min_step)

            # Create levels
            levels_hw = np.arange(min_val, max_val + rounded_step, rounded_step)
            n_levels = len(levels_hw)

            # Clip if too many levels
            if n_levels > max_n_levels:
                rounded_step = range_val / max_n_levels
                levels_hw = np.arange(min_val, max_val + rounded_step, rounded_step)

            # Recompute levels if not enough
            if len(levels_hw) < 3:
                levels_hw = np.linspace(min_val, max_val, num=3)

            # Create ticks aligned with levels
            tick_step = range_val / target_n_ticks
            rounded_tick_step = max(np.round(tick_step, 2), min_step)

            ticks_z = np.arange(min_val, max_val + rounded_tick_step, rounded_tick_step)
            ticks_z = np.round(ticks_z, 2)

            # Clip ticks if too many
            if len(ticks_z) > max_n_ticks:
                ticks_z = np.linspace(min_val, max_val, num=max_n_ticks)
                ticks_z = np.round(ticks_z, 2)

            # Define colorbar settings
            cbar_kwargs = {'ticks': ticks_z}

            # Plot the filled contour
            p = da.plot.contourf(
                x='birth_year',
                y='GMT',
                ax=ax,
                add_labels=False,
                levels=levels_hw,
                cmap='Reds',
                extend='both',
                cbar_kwargs=cbar_kwargs
            )

        ax.set_yticks(
            ticks=GMT_indices_ticks,
            labels=gmts2100,
            color='black',
        )

        ax.set_xticks(
            ticks=np.arange(1960, 2025, 10),
            color='black',
        )

        ax.tick_params(axis='both', labelsize=14)

        # Title and aesthetics
        ax.set_title(
        "Lifetime Exposure to {} \n in {}".format(extr_name, country),
        fontsize=18, fontweight='bold', loc='center', color='black'
        )

        ax.spines['right'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')

        ax.set_xlabel('Birth year', fontsize=16, color='black')
        ax.set_ylabel('GMT warming by 2100 [°C]', fontsize=16, color='black', rotation='vertical')

        plt.tight_layout()

        if flags['old_demo']:

            plt.savefig(
                scripts_dir + '/figures/source2suffering/sandbox/old_demography/fig8_BE_{}_lifetime_exposure_{}_gmt_{}_{}.png'.format(
                    extr, country, flags['gmt'],flags['rm']
                )
            )

        else:

            plt.savefig(
                scripts_dir + '/figures/source2suffering/development/{}/fig8_BE_{}_lifetime_exposure_{}_gmt_{}_{}.png'.format(
                    extr, extr, country, flags['gmt'],flags['rm']
                )
            )

    print(f"Plot performed for {country}")


def plot_dev_fig9(
    s2s_valc_nr_children_facing_extra_hazard,
    wt_valc_nr_children_facing_extra_hazard,
    birth_cohort_int,
    hazards,
    flags,
):
    """
    Plot comparison of the assessment report produced by W.Thiery et al.(2021) and the Source2Suffering framework 

    Parameters:
    - s2s_valc_nr_children_facing_extra_hazard (Array) : Contains per birth cohort of interested the number of extra children affected computed by the Source2Suffering framework
    - wt_valc_nr_children_facing_extra_hazard (Array) : Contains per birth cohort of interested the number of extra children affected computed by W.Thiery
    - birth_cohort_int (Array) : Contains the birth cohort of interest
    - flags (Array) : contains the parameters of the configuration
    """

    plt.close('all')

    plt.figure(figsize=(13, 9))

    if hazards=='burntarea':
        hazards_name = 'Wildfires'
    if hazards=='cropfailedarea':
        hazards_name = 'Crop failures'
    if hazards=='driedarea':
        hazards_name = 'Droughts'
    if hazards=='floodedarea':
        hazards_name = 'River floods'
    if hazards=='heatwavedarea':
        hazards_name = 'Heatwaves'
    if hazards=='tropicalcyclonedarea':
        hazards_name = 'Tropical Cyclones'

    # Plot lines with markers (instead of only points)
    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard, 
            '-o', color="tab:orange", label="S2S flags['gmt'] = {} & flags['rm'] = {}".format(flags['gmt'],flags['rm']), markersize=10,lw=2)
    plt.plot(birth_cohort_int, wt_valc_nr_children_facing_extra_hazard, 
            '-o', color="tab:blue", label="W.Thiery", markersize=10,lw=2)

    # Set the figure title
    plt.title("Number of children facing at\nleast one extra {} due to NeptunDeep".format(hazards_name),
            fontsize=22, fontweight='bold')

    # Axis labels
    plt.xlabel("Birth cohort", fontsize=20)
    plt.ylabel("Number of children", fontsize=20)

    # Display the legend
    plt.legend(fontsize=18)

    # Set x-axis ticks for each birth cohort
    plt.xticks(birth_cohort_int, rotation=45, fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid()

    plt.tight_layout()

    # Save the figure
    plt.savefig(scripts_dir + '/figures/assessment/NeptunDeep/fig9_nr_children_facing_extra_{}_NeptunDeep_gmt_{}_{}.png'.format(hazards, flags['gmt'],flags['rm']))

def plot_dev_fig10(
    s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new,
    s2s_valc_nr_children_facing_extra_hazard_gmt_original,
    wt_valc_nr_children_facing_extra_hazard,
    birth_cohort_int,
    hazards,
    flags,
):
    """
    Plot comparison of the assessment report produced by W.Thiery et al.(2021) and the Source2Suffering framework and the two different GMT options

    Parameters:
    - s2s_valc_nr_children_facing_extra_hazard (Array) : Contains per birth cohort of interested the number of extra children affected computed by the Source2Suffering framework
    - wt_valc_nr_children_facing_extra_hazard (Array) : Contains per birth cohort of interested the number of extra children affected computed by W.Thiery
    - birth_cohort_int (Array) : Contains the birth cohort of interest
    - flags (Array) : contains the parameters of the configuration
    """

    plt.close('all')

    plt.figure(figsize=(13, 9))

    if hazards=='burntarea':
        hazards_name = 'Wildfires'
    if hazards=='cropfailedarea':
        hazards_name = 'Crop failures'
    if hazards=='driedarea':
        hazards_name = 'Droughts'
    if hazards=='floodedarea':
        hazards_name = 'River floods'
    if hazards=='heatwavedarea':
        hazards_name = 'Heatwaves'
    if hazards=='tropicalcyclonedarea':
        hazards_name = 'Tropical Cyclones'

    # Plot lines with markers (instead of only points)

    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new, 
            '-o', color="tab:orange", label="S2S - GMT AR6_new & flags['rm'] = {}".format(flags['rm']), markersize=10,lw=2)

    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard_gmt_original, 
            '-o', color="tab:red", label="S2S - GMT Original & flags['rm'] = {}".format(flags['rm']), markersize=10,lw=2)

    plt.plot(birth_cohort_int, wt_valc_nr_children_facing_extra_hazard, 
            '-o', color="tab:blue", label="W.Thiery", markersize=10,lw=2)

    # Set the figure title
    plt.title("Number of children facing at\nleast one extra {} due to NeptunDeep".format(hazards_name),
            fontsize=22, fontweight='bold')

    # Axis labels
    plt.xlabel("Birth cohort", fontsize=20)
    plt.ylabel("Number of children", fontsize=20)

    # Display the legend
    plt.legend(fontsize=18)

    # Set x-axis ticks for each birth cohort
    plt.xticks(birth_cohort_int, rotation=45, fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid()

    plt.tight_layout()

    # Save the figure
    plt.savefig(scripts_dir + '/figures/assessment/NeptunDeep/fig10_nr_children_facing_extra_{}_NeptunDeep_gmt_comp_{}.png'.format(hazards,flags['rm']))


def plot_dev_fig11(df_GMT_15, df_GMT_20, df_GMT_NDC, df_GMT_OS, df_GMT_noOS):
    """
    Plot the time evolution of global mean temperature trajectories for different scenarios,
    limited to years up to 2100.
    """
    plt.figure(figsize=(10, 6))

    # Plot each trajectory up to year 2100
    plt.plot(df_GMT_15.loc[:2100].index, df_GMT_15.loc[:2100].values, label='1.5°C', color='tab:blue', lw=3)
    plt.plot(df_GMT_20.loc[:2100].index, df_GMT_20.loc[:2100].values, label='2.0°C', color='tab:green', lw=3)
    plt.plot(df_GMT_NDC.loc[:2100].index, df_GMT_NDC.loc[:2100].values, label='NDC', color='tab:orange', lw=3)
    plt.plot(df_GMT_OS.loc[:2100].index, df_GMT_OS.loc[:2100].values, label='OS', color='tab:red', lw=3)
    plt.plot(df_GMT_noOS.loc[:2100].index, df_GMT_noOS.loc[:2100].values, label='noOS', color='tab:purple', lw=3)

    # Add axis labels and title
    plt.xlabel('Year',fontsize=20)
    plt.ylabel('GMT anomaly (°C)',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Global Mean Temperature Trajectories (up to 2100)',fontsize=22)
    plt.legend(fontsize=18)
    plt.grid(True)
    plt.tight_layout()

    # Save the figure
    plt.savefig(scripts_dir + '/figures/source2suffering/development/GMT_trajectories/GMT_traj_thiery_et_al.png')

def plot_dev_fig12(df_GMT_OS, df_GMT_noOS, ds_GMT_STS):
    """
    Plot the time evolution of global mean temperature trajectories for different scenarios,
    limited to years up to 2100.
    """
    plt.figure(figsize=(10, 6))

    # unpack the values of the four differents GMT pathways inside the STS scenarios #

    da_GMT_STS_ModAct = ds_GMT_STS['tas'].sel(
    time=slice(1960, 2100),
    percentile='50.0',
    scenario='ModAct'
    )

    da_GMT_STS_Ren = ds_GMT_STS['tas'].sel(
    time=slice(1960, 2100),
    percentile='50.0',
    scenario='Ren'
    )

    da_GMT_STS_LD = ds_GMT_STS['tas'].sel(
    time=slice(1960, 2100),
    percentile='50.0',
    scenario='LD'
    )

    da_GMT_STS_SP = ds_GMT_STS['tas'].sel(
    time=slice(1960, 2100),
    percentile='50.0',
    scenario='SP'
    )

    # Plot each trajectory up to year 2100
    plt.plot(df_GMT_OS.loc[:2100].index, df_GMT_OS.loc[:2100].values, label='OS', color='tab:red', lw=3)
    plt.plot(df_GMT_noOS.loc[:2100].index, df_GMT_noOS.loc[:2100].values, label='noOS', color='tab:purple', lw=3)
    plt.plot(df_GMT_noOS.loc[:2100].index, da_GMT_STS_ModAct, label='ModAct', color='tab:blue', lw=3, linestyle="--")
    plt.plot(df_GMT_noOS.loc[:2100].index, da_GMT_STS_Ren, label='Ren', color='tab:cyan', lw=3, linestyle="--")

    # Add axis labels and title
    plt.xlabel('Year',fontsize=20)
    plt.ylabel('GMT anomaly (°C)',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Global Mean Temperature Trajectories (up to 2100)',fontsize=22)
    plt.legend(fontsize=18)
    plt.grid(True)
    plt.tight_layout()

    # Save the figure
    plt.savefig(scripts_dir + '/figures/source2suffering/development/GMT_trajectories/GMT_traj_STS_Thiery_comp.png')

def plot_dev_fig13_regions(ds_regions, extr, flags, ds_le, region_ind, EMF):
    """
    Plot MMM lifetime exposure to heatwave for all regions, and the two STS scenarios for the SPARCCLE Delivery.

    Parameters:
    - ds_regions (xarray.Dataset): Dataset with all the countries per region
    - ds_le (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - region (int): Index of the region (ex: '1')
    """
    
    plt.close('all') 

    plt.figure(figsize=(12, 8))

    GMT_color = ['tab:red','tab:blue']
    GMT_label = ['ModAct', 'Ren']  
    region_name = ds_regions['name'].sel(region=region_ind)

    if extr=='burntarea':
        extr_name = 'Wildfires'
    if extr=='cropfailedarea':
        extr_name = 'Crop failures'
    if extr=='driedarea':
        extr_name = 'Droughts'
    if extr=='floodedarea':
        extr_name = 'River floods'
    if extr=='heatwavedarea':
        extr_name = 'Heatwaves'
    if extr=='tropicalcyclonedarea':
        extr_name = 'Tropical Cyclones'

    le_ModAct = ds_le['mmm_STS_ModAct'].sel(
        region=region_ind,
    )

    le_ModAct_std = ds_le['std_STS_ModAct'].sel(
        region=region_ind,
    )/2
    
    le_Ren = ds_le['mmm_STS_Ren'].sel(
        region=region_ind,
    )

    le_Ren_std = ds_le['std_STS_Ren'].sel(
        region=region_ind,
    )/2

    plt.plot(birth_years, le_ModAct, linestyle='-',color=GMT_color[0],label=GMT_label[0],lw=3)
    plt.fill_between(birth_years, le_ModAct - le_ModAct_std, le_ModAct + le_ModAct_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[0], lw=3)
    plt.plot(birth_years, le_Ren, linestyle='-',color=GMT_color[1],label=GMT_label[1],lw=3)
    plt.fill_between(birth_years, le_Ren - le_Ren_std, le_Ren + le_Ren_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[1], lw=3)

    plt.title("Lifetime Exposure to {} \n in the {} region".format(extr_name, region_name),fontsize=18,fontweight='bold')
    plt.xlabel("Birth Year",fontsize=16)
    plt.ylabel("Lifetime Exposure",fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=16)
    plt.savefig(scripts_dir+'/figures/assessment/SPARCCLE_STS/lifetime_exposure_{}_region_{}_{}.png'.format(extr,region_name,flags['rm']))

def plot_dev_fig14_regions(flags, extr, ds_lfe_perregion, ind_region):
    """
    Plot MMM fraction of land exposed to hazard for the a specific region under the Thiery et al.(2021) trajectories.

    Parameters:
    - extr (str): flags['extr']
    - ds_lfe_perregion (xarray.Dataset): Dataset with all the countries per region
    - ind_region (int): Index of the region under interest
    """
    
    plt.close('all') 

    plt.figure(figsize=(12, 8))

    if extr=='burntarea':
        extr_name = 'Wildfires'
    if extr=='cropfailedarea':
        extr_name = 'Crop failures'
    if extr=='driedarea':
        extr_name = 'Droughts'
    if extr=='floodedarea':
        extr_name = 'River floods'
    if extr=='heatwavedarea':
        extr_name = 'Heatwaves'
    if extr=='tropicalcyclonedarea':
        extr_name = 'Tropical Cyclones'

    if flags['gmt']=='original':

        GMT_label = ['1.5°C', '2.0°C','NDC']  

    if flags['gmt']=='ar6_new':

        GMT_list = [0,11,20]
        GMT_label = ['1.5°C', '2.5°C','3.5°C']  

    GMT_color = ['#FFC72C', '#FF5A1F', '#C74135']
    region_name = ds_regions['name'].sel(region=ind_region)

    if flags['gmt']=='original':

        lfe_15 = ds_lfe_perregion['mmm_15_sm'].sel(region=ind_region).isel(time_ind=slice(0, -13))

        y_std_15 = ds_lfe_perregion['std_15_sm'].sel(region=ind_region).isel(time_ind=slice(0, -13))/2

        lfe_20 = ds_lfe_perregion['mmm_20_sm'].sel(region=ind_region).isel(time_ind=slice(0, -13))

        y_std_20 = ds_lfe_perregion['std_20_sm'].sel(region=ind_region).isel(time_ind=slice(0, -13))/2

        lfe_NDC = ds_lfe_perregion['mmm_NDC_sm'].sel(region=ind_region).isel(time_ind=slice(0, -13))

        y_std_NDC = ds_lfe_perregion['std_NDC_sm'].sel(region=ind_region).isel(time_ind=slice(0, -13))/2

        plt.plot(year_range[:-13], lfe_15*100, linestyle='-',color=GMT_color[0],label=GMT_label[0], lw=3)
        plt.fill_between(year_range[:-13], lfe_15*100 - y_std_15*100, lfe_15*100 + y_std_15*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[0], lw=3)

        plt.plot(year_range[:-13], lfe_20*100, linestyle='-',color=GMT_color[1],label=GMT_label[1], lw=3)
        plt.fill_between(year_range[:-13], lfe_20*100 - y_std_20*100, lfe_20*100 + y_std_20*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[1], lw=3)

        plt.plot(year_range[:-13],lfe_NDC*100, linestyle='-',color=GMT_color[2],label=GMT_label[2], lw=3)
        plt.fill_between(year_range[:-13], lfe_NDC*100 - y_std_NDC*100, lfe_NDC*100 + y_std_NDC*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[2], lw=3)

    # Plot script works but I have not yet generelized the pipeline before to produce the object of the mmm_BE_sm #

    if flags['gmt']=='ar6_new':
        
        print("\nPlot not performed for the {} region because routine with the stylized trajectories not yet settle for the LFE with flags['gmt'] = ar6_new".format(region_name))

        return 

    #         for GMT in GMT_list:

    #             lfe = ds_lfe_perregion['mmm_BE_sm'].sel(
    #                 region=ind_region,
    #                 GMT=GMT
    #             ).isel(time_ind=slice(0, -13))

    #             y_std = ds_lfe_perregion['std_BE_sm'].sel(
    #                 region=ind_region,
    #                 GMT = GMT
    #             ).isel(time_ind=slice(0, -13))/2

    #             plt.plot(year_range[:-13], lfe*100, linestyle='-',color=GMT_color[i],label=GMT_label[i],lw=3)
    #             plt.fill_between(year_range[:-13], lfe*100 - y_std*100, lfe*100 + y_std*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i],lw=3)

    plt.title("Land Fraction Exposed to {} \n in the {} region".format(extr_name, region_name),fontsize=18,fontweight='bold')
    plt.xlabel("Time",fontsize=16,fontweight='bold')
    plt.ylabel("Area annually exposed to {} (%)".format(extr_name),fontsize=16,fontweight='bold')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=16)
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig14_mmm_landfrac_exposed_region_{}_gmt_{}_{}.png'.format(extr,region_name,flags['gmt'],flags['rm']))
    print("\nPlot performed for the {} region".format(region_name))

def plot_dev_fig15_country(flags, extr, ds_lfe_percountry, country_name):
    """
    Plot MMM fraction of land exposed to hazard for all countries under the Thiery et al.(2021) trajectories.

    Parameters:
    - extr (str): flags['extr']
    - ds_lfe_percountry (xarray.Dataset): Dataset with all the countries 
    - country (str): Name of the country under study 
    """
    
    plt.close('all') 

    plt.figure(figsize=(12, 8))

    if extr=='burntarea':
        extr_name = 'Wildfires'
    if extr=='cropfailedarea':
        extr_name = 'Crop failures'
    if extr=='driedarea':
        extr_name = 'Droughts'
    if extr=='floodedarea':
        extr_name = 'River floods'
    if extr=='heatwavedarea':
        extr_name = 'Heatwaves'
    if extr=='tropicalcyclonedarea':
        extr_name = 'Tropical Cyclones'

    if flags['gmt']=='original':

        GMT_label = ['1.5°C', '2.0°C','NDC']  

    if flags['gmt']=='ar6_new':

        GMT_list = [0,11,20]
        GMT_label = ['1.5°C', '2.5°C','3.5°C']  

    GMT_color = ['#FFC72C', '#FF5A1F', '#C74135']

    if flags['gmt']=='original':

        lfe_15 = ds_lfe_percountry['mmm_15_sm'].sel(country=country_name).isel(time_ind=slice(0, -13))

        y_std_15 = ds_lfe_percountry['std_15_sm'].sel(country=country_name).isel(time_ind=slice(0, -13))/2

        lfe_20 = ds_lfe_percountry['mmm_20_sm'].sel(country=country_name).isel(time_ind=slice(0, -13))

        y_std_20 = ds_lfe_percountry['std_20_sm'].sel(country=country_name).isel(time_ind=slice(0, -13))/2

        lfe_NDC = ds_lfe_percountry['mmm_NDC_sm'].sel(country=country_name).isel(time_ind=slice(0, -13))

        y_std_NDC = ds_lfe_percountry['std_NDC_sm'].sel(country=country_name).isel(time_ind=slice(0, -13))/2

        plt.plot(year_range[:-13], lfe_15*100, linestyle='-',color=GMT_color[0],label=GMT_label[0], lw=3)
        plt.fill_between(year_range[:-13], lfe_15*100 - y_std_15*100, lfe_15*100 + y_std_15*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[0], lw=3)

        plt.plot(year_range[:-13], lfe_20*100, linestyle='-',color=GMT_color[1],label=GMT_label[1], lw=3)
        plt.fill_between(year_range[:-13], lfe_20*100 - y_std_20*100, lfe_20*100 + y_std_20*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[1], lw=3)

        plt.plot(year_range[:-13],lfe_NDC*100, linestyle='-',color=GMT_color[2],label=GMT_label[2], lw=3)
        plt.fill_between(year_range[:-13], lfe_NDC*100 - y_std_NDC*100, lfe_NDC*100 + y_std_NDC*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[2], lw=3)

    # Plot script works but I have not yet generelized the pipeline before to produce the object of the mmm_BE_sm #

    if flags['gmt']=='ar6_new':
        
        #print("\nPlot not performed for {} because routine with the stylized trajectories not yet settle for the LFE with flags['gmt'] = ar6_new".format(country_name))
        
        #return
            
        for GMT in GMT_list:

            lfe = ds_lfe_percountry['mmm_BE_sm'].sel(
                country=country_name,
                GMT=GMT
            ).isel(time_ind=slice(0, -13))

            y_std = ds_lfe_percountry['std_BE_sm'].sel(
                country=country_name,
                GMT = GMT
            ).isel(time_ind=slice(0, -13))

            plt.plot(year_range[:-13], lfe*100, linestyle='-',color=GMT_color[i],label=GMT_label[i],lw=3)
            plt.fill_between(year_range[:-13], lfe*100 - y_std*100, lfe*100 + y_std*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i],lw=3)

    plt.title("Land Fraction Exposed to {} \n in {}".format(extr_name, country_name),fontsize=18,fontweight='bold')
    plt.xlabel("Time",fontsize=16,fontweight='bold')
    plt.ylabel("Area annually exposed to {} (%)".format(extr_name),fontsize=16,fontweight='bold')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=16)
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig15_mmm_landfrac_exposed_{}_gmt_{}_{}.png'.format(extr,country_name,flags['gmt'],flags['rm']))
    print("Plot performed for {}".format(country_name))

def plot_dev_fig16(
    wt_valc_nr_children_facing_extra_hazard,
    s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new_rm,
    s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new_no_rm,
    s2s_valc_nr_children_facing_extra_hazard_gmt_original_rm,
    s2s_valc_nr_children_facing_extra_hazard_gmt_original_no_rm,
    birth_cohort_int,
    hazards,
):
    """
    Plot comparison of the assessment report produced by W.Thiery et al.(2021) and the Source2Suffering framework and the two different GMT options + the thwo different rm options

    Parameters:
    - wt_valc_nr_children_facing_extra_hazard (Array) : Contains per birth cohort of interested the number of extra children affected computed by W.Thiery
    - s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new_rm (Array) : Value produce by the S2S framework with the ar6_new and rm configuration
    - s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new_no_rm (Array) : Value produce by the S2S framework with the ar6_new and no_rm configuration
    - s2s_valc_nr_children_facing_extra_hazard_gmt_original_rm (Array) : Value produce by the S2S framework with the original and rm configuration
    - s2s_valc_nr_children_facing_extra_hazard_gmt_original_no_rm (Array) : Value produce by the S2S framework with the original and no_rm configuration
    - birth_cohort_int (Array) : Contains the birth cohort of interest
    - flags (Array) : contains the parameters of the configuration
    """

    plt.close('all')

    plt.figure(figsize=(13, 9))

    if hazards=='burntarea':
        hazards_name = 'Wildfires'
    if hazards=='cropfailedarea':
        hazards_name = 'Crop failures'
    if hazards=='driedarea':
        hazards_name = 'Droughts'
    if hazards=='floodedarea':
        hazards_name = 'River floods'
    if hazards=='heatwavedarea':
        hazards_name = 'Heatwaves'
    if hazards=='tropicalcyclonedarea':
        hazards_name = 'Tropical Cyclones'

    # Plot lines with markers (instead of only points)


    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard_gmt_original_rm, 
            '--o', color="tab:orange", label="S2S - GMT Original & {}".format('rm'), markersize=10,lw=3)

    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new_no_rm, 
            '-o', color="tab:red", label="S2S - GMT AR6_new & {}".format('no_rm'), markersize=10,lw=3)

    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard_gmt_ar6_new_rm, 
            '--o', color="tab:red", label="S2S - GMT AR6_new & {}".format('rm'), markersize=10,lw=3)

    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard_gmt_original_no_rm, 
            '-o', color="tab:orange", label="S2S - GMT Original & {}".format('no_rm'), markersize=10,lw=6)

    plt.plot(birth_cohort_int, wt_valc_nr_children_facing_extra_hazard, 
            '-o', color="tab:blue", label="W.Thiery", markersize=10,lw=3)

    # Set the figure title
    plt.title("Number of children facing at\nleast one extra {} due to NeptunDeep".format(hazards_name),
            fontsize=22, fontweight='bold')

    # Axis labels
    plt.xlabel("Birth cohort", fontsize=20, fontweight='bold')
    plt.ylabel("Number of children", fontsize=20, fontweight='bold')

    # Display the legend
    plt.legend(fontsize=18)

    # Set x-axis ticks for each birth cohort
    plt.xticks(birth_cohort_int, rotation=45, fontsize=14)
    plt.yticks(fontsize=14, rotation=45)
    plt.grid()

    plt.tight_layout()

    # Save the figure
    plt.savefig(scripts_dir + '/figures/assessment/NeptunDeep/fig16_nr_children_facing_extra_{}_NeptunDeep_gmt_and_rm_comp.png'.format(hazards))

def plot_dev_fig17_regions(ds_regions, extr, flags, ds_lfe, region_ind, EMF):
    """
    Plot MMM landfraction exposed to specific hazard for all regions, and the two STS scenarios for the SPARCCLE Delivery.

    Parameters:
    - ds_regions (xarray.Dataset): Dataset with all the countries per region
    - ds_lfe (xarray.Dataset): Dataset containing the "land fraction exposed" variable
    - region (int): Index of the region (ex: '1')
    """
    
    plt.close('all') 

    plt.figure(figsize=(12, 8))

    GMT_color = ['tab:red','tab:blue']
    GMT_label = ['ModAct', 'Ren']  
    region_name = ds_regions['name'].sel(region=region_ind)

    if extr=='burntarea':
        extr_name = 'Wildfires'
    if extr=='cropfailedarea':
        extr_name = 'Crop failures'
    if extr=='driedarea':
        extr_name = 'Droughts'
    if extr=='floodedarea':
        extr_name = 'River floods'
    if extr=='heatwavedarea':
        extr_name = 'Heatwaves'
    if extr=='tropicalcyclonedarea':
        extr_name = 'Tropical Cyclones'

    lfe_ModAct = ds_lfe['mmm_STS_ModAct_sm'].sel(
        region=region_ind,
    )

    lfe_ModAct_std = ds_lfe['std_STS_ModAct_sm'].sel(
        region=region_ind,
    )/2
    
    lfe_Ren = ds_lfe['mmm_STS_Ren_sm'].sel(
        region=region_ind,
    )

    lfe_Ren_std = ds_lfe['std_STS_Ren_sm'].sel(
        region=region_ind,
    )/2

    plt.plot(year_range[:-13], lfe_ModAct.isel(time_ind=slice(0, -13))*100, linestyle='-',color=GMT_color[0],label=GMT_label[0],lw=3)
    plt.fill_between(year_range[:-13], (lfe_ModAct.isel(time_ind=slice(0, -13)) - lfe_ModAct_std.isel(time_ind=slice(0, -13)))*100, (lfe_ModAct.isel(time_ind=slice(0, -13)) + lfe_ModAct_std.isel(time_ind=slice(0, -13)))*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[0], lw=3)
    plt.plot(year_range[:-13], lfe_Ren.isel(time_ind=slice(0, -13))*100, linestyle='-',color=GMT_color[1],label=GMT_label[1],lw=3)
    plt.fill_between(year_range[:-13], (lfe_Ren.isel(time_ind=slice(0, -13)) - lfe_Ren_std.isel(time_ind=slice(0, -13)))*100, (lfe_Ren.isel(time_ind=slice(0, -13)) + lfe_Ren_std.isel(time_ind=slice(0, -13)))*100, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[1], lw=3)

    plt.title("Land Fraction Exposed to {} \n in {}".format(extr_name, region_name),fontsize=18,fontweight='bold')
    plt.xlabel("Time",fontsize=16,fontweight='bold')
    plt.ylabel("Area annually exposed to {} (%)".format(extr_name),fontsize=16,fontweight='bold')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=16)
    plt.savefig(scripts_dir+'/figures/assessment/SPARCCLE_STS/landfraction_exposed_{}_region_{}.png'.format(extr,region_name))