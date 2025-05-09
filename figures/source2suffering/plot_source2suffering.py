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

def plot_dev_fig1(flags, ds, country, run, GMT):
    """
    Plot lifetime exposure for a specific country, run, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - run (int): Numéro de simulation (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """

    plt.close('all') 
    
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
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig1_lifetime_exposure_{}_GMT_{}_isimip_sim{}.png'.format(flags['extr'],country,GMT,run))

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
        filename = scripts_dir + '/figures/source2suffering/development/{}/fig2_lifetime_exposure_trends_isimip_sim{}_GMT_{}_region_{}.png'.format(flags['extr'],
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
    
    plt.close('all') 

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
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig3_lifetime_exposure_region_{}_GMT_{}_isimip_sim{}.png'.format(flags['extr'],region_name,GMT,run))


def plot_dev_fig4(flags, ds, country, GMT):
    """
    Plot MMM lifetime exposure for a specific country, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (str): Nom du pays (ex: 'Belgium')
    - GMT (str or int): Niveau GMT (ex: '21')
    """

    plt.close('all') 
    
    exposure = ds['mmm'].sel(
        country=country,
        GMT=GMT  
    )

    y_std = ds['std'].sel(
        country = country,
        GMT = GMT
    )

    birth_years = ds['birth_year'].values

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$')
    plt.title("{} : MMM Lifetime Exposure in {} (GMT {})".format(flags['extr'], country, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig4_mmm_lifetime_exposure_{}_GMT_{}.png'.format(flags['extr'],country,GMT))


def plot_dev_fig5(ds_regions, flags, ds, region, GMT):
    """
    Plot MMM lifetime exposure for a specific region, and GMT value.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - region (int): Region index (ex: 6)
    - GMT (str or int): Niveau GMT (ex: '21')
    """
    
    plt.close('all') 

    exposure = ds['mmm'].sel(
        region = region,
        GMT = GMT  
    )

    birth_years = ds['birth_year'].values
    region_name = ds_regions['name'].sel(region=region)
    y_std = ds['std'].sel(
        region = region,
        GMT = GMT
    )

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(birth_years, exposure, marker='o', linestyle='-')
    plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$')
    plt.title("{} : MMM Lifetime Exposure in {} (GMT {})".format(flags['extr'], region_name, GMT))
    plt.xlabel("Birth Year")
    plt.ylabel("Lifetime Exposure")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig5_mmm_lifetime_exposure_region_{}_GMT_{}.png'.format(flags['extr'],region_name,GMT))

def plot_dev_fig6_country(flags, ds, country, EMF):
    """
    Plot MMM lifetime exposure to heatwave for all countries, and three GMT values.

    Parameters:
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - country (string): Name of the country (ex: 'Belgium')
    - GMT (str or int): Niveau GMT (ex: '21')
    - EMF (Boolean) : Chose to express the results in terms of the EMF or remains in the absolute lifetime exposure
    """
    
    plt.close('all') 

    plt.figure(figsize=(10, 6))

    GMT_list = [0,11,20]

    GMT_color = ['tab:blue','tab:orange','tab:red']
    GMT_label = ['1.5°C', '2.5°C', '3.5°C']  
    i = 0

    if EMF:

        for GMT in GMT_list:
    
            exposure = ds.sel(
                country=country,
                GMT=GMT  
            )

            # y_std = ds['std'].sel(
            #     country = country,
            #     GMT = GMT
            # )

            plt.plot(birth_years, exposure, linestyle='-',color=GMT_color[i],label=GMT_label[i])
            #plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i])

            i +=1

        plt.title("{} : MMM Lifetime Exposure EMF in {}".format(flags['extr'], country))
        plt.xlabel("Birth Year")
        plt.ylabel("Exposure Multiplication Factor (EMF)")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/storyline/fig6_mmm_lifetime_exposure_EMF_{}.png'.format(flags['extr'],country))

    else:

        for GMT in GMT_list:

            exposure = ds['mmm'].sel(
                country=country,
                GMT=GMT  
            )

            y_std = ds['std'].sel(
                country = country,
                GMT = GMT
            )

            plt.plot(birth_years, exposure, linestyle='-',color=GMT_color[i],label=GMT_label[i])
            plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i])

            i +=1

        plt.title("{} : MMM Lifetime Exposure in {}".format(flags['extr'], country))
        plt.xlabel("Birth Year")
        plt.ylabel("Lifetime Exposure")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/storyline/fig6_mmm_lifetime_exposure_{}.png'.format(flags['extr'],country))

def plot_dev_fig6_region(ds_regions, flags, ds, region, EMF):
    """
    Plot MMM lifetime exposure to heatwave for all regions, and three GMT values.

    Parameters:
    - ds_regions (xarray.Dataset): Dataset with all the countries per region
    - ds (xarray.Dataset): Dataset contenant la variable 'lifetime_exposure'
    - region (int): Index of the region (ex: '1')
    """
    
    plt.close('all') 

    plt.figure(figsize=(10, 6))

    GMT_list = [0,11,20]

    GMT_color = ['tab:blue','tab:orange','tab:red']
    GMT_label = ['1.5°C', '2.5°C', '3.5°C']  
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

        plt.title("{} : MMM Lifetime Exposure EMF in {}".format(flags['extr'], region_name))
        plt.xlabel("Birth Year")
        plt.ylabel("Exposure Multiplication Factor (EMF)")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/storyline/fig6_mmm_lifetime_exposure_EMF_region_{}.png'.format(flags['extr'],region_name))

    else: 

        for GMT in GMT_list:

            exposure = ds['mmm'].sel(
                region=region,
                GMT=GMT  
            )

            y_std = ds['std'].sel(
                region=region,
                GMT = GMT
            )

            plt.plot(birth_years, exposure, linestyle='-',color=GMT_color[i],label=GMT_label[i])
            plt.fill_between(birth_years, exposure - y_std, exposure + y_std, alpha=0.3, label=r'$\pm \sigma$',color=GMT_color[i])

            i +=1

        plt.title("{} : MMM Lifetime Exposure in {}".format(flags['extr'], region_name))
        plt.xlabel("Birth Year")
        plt.ylabel("Lifetime Exposure")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/storyline/fig6_mmm_lifetime_exposure_region_{}.png'.format(flags['extr'],region_name))

def plot_dev_fig7(
    ds_regions,
    df_GMT_strj,
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
    GMT_indices_ticks=[0,5,10,15,20]
    gmts2100 = np.round(df_GMT_strj.loc[2100,GMT_indices_ticks].values,1)    
 
    
    #--------------------------- Load Data --------------------------------#

    with open(data_dir+'{}/{}/lifetime_exposure_perregion_GMT.pkl'.format(flags['version'],flags['extr']), 'rb') as f:
            ds_le_perregion_GMT = pk.load(f)

    #------------------------------ Plot ----------------------------------#

    # Ticks style

    mpl.rcParams['xtick.labelcolor'] = 'gray'
    mpl.rcParams['ytick.labelcolor'] = 'gray'

    # Size of the figure

    fig = plt.figure(figsize=(8, 7))

    # Init one axe 

    ax = plt.gca()

    # Main plot

    da = ds_le_perregion_GMT['mmm'].loc[{
        'region':region_ind,
        'birth_year': np.arange(1960, 2021),
    }]


    if flags['extr']=='heatwavedarea':

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

        # Compute the maximum value of the data
        max_val = da.max(skipna=True).item()

        # Set a target and maximum number of contour levels
        target_n_levels = 30
        max_n_levels = 50

        # Compute raw step
        step = max_val / target_n_levels if target_n_levels != 0 else 0.1

        # Add a minimum threshold for the step to avoid too many levels
        min_step = 0.1
        rounded_step = max(np.round(step, 2), min_step)

        # Compute levels, but clip number of levels to avoid overflow
        n_levels = int(np.ceil(max_val / rounded_step)) + 1

        # Safety check to avoid exceeding the max allowed size
        if n_levels > max_n_levels:
            rounded_step = max_val / max_n_levels
            n_levels = max_n_levels

        # Create levels array
        levels_hw = np.arange(0, rounded_step * n_levels, rounded_step)

        # Target and maximum number of ticks for the colorbar (for readability)
        target_n_ticks = 6
        max_n_ticks = 10

        # Compute a raw tick step
        tick_step = max_val / target_n_ticks if target_n_ticks != 0 else 0.5

        # Add a minimum threshold for the tick step
        min_tick_step = 0.1
        rounded_tick_step = max(np.round(tick_step, 2), min_tick_step)

        # Compute number of ticks, enforce a max limit
        n_ticks = int(np.ceil(max_val / rounded_tick_step)) + 1

        if n_ticks > max_n_ticks:
            rounded_tick_step = max_val / max_n_ticks
            n_ticks = max_n_ticks

        # Generate tick positions safely
        ticks_z = np.arange(0, rounded_tick_step * n_ticks, rounded_tick_step)

        ticks_z = np.round(ticks_z, 1)

        # Use in colorbar kwargs
        cbar_kwargs = {'ticks': ticks_z}

        p =da.plot.contourf(
            x='birth_year',
            y='GMT',
            ax=ax,
            add_labels=False,
            levels=levels_hw,
            cmap='Reds',
            cbar_kwargs={'ticks': ticks_z}
        )

    ax.set_yticks(
            ticks=GMT_indices_ticks,
            labels=gmts2100,
            color='gray',
    )
    
    ax.set_xticks(
        ticks=np.arange(1960,2025,10),
        color='gray',
    )  

    ax.tick_params(axis='both', labelsize=12)  

    ax.set_title(
        "{} : Lifetime Exposure in \n {}".format(flags['extr'], region_name),
        loc='center',
        fontweight='bold',
        color='gray',
        fontsize=16,
    )
    ax.spines['right'].set_color('gray')
    ax.spines['top'].set_color('gray')
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_color('gray')  
        
    ax.set_xlabel('Birth year',fontsize=14,color='gray')
    ax.set_ylabel('GMT warming by 2100 [°C]',fontsize=14,color='gray',rotation='vertical',)

    plt.tight_layout()
    
    plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig7_BE_{}_lifetime_exposure_region_{}.png'.format(flags['extr'],flags['extr'],region_name))


def plot_dev_fig8(
    df_GMT_strj,
    flags,
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
    GMT_indices_ticks=[0,5,10,15,20]
    gmts2100 = np.round(df_GMT_strj.loc[2100,GMT_indices_ticks].values,1)    
 
    
    #--------------------------- Load Data --------------------------------#

    with open(data_dir+'{}/{}/lifetime_exposure_percountry_GMT.pkl'.format(flags['version'],flags['extr']), 'rb') as f:
            ds_le_percountry_GMT = pk.load(f)

    #------------------------------ Plot ----------------------------------#

    # Ticks style

    mpl.rcParams['xtick.labelcolor'] = 'gray'
    mpl.rcParams['ytick.labelcolor'] = 'gray'

    # Size of the figure

    fig = plt.figure(figsize=(8, 7))

    # Init one axe 

    ax = plt.gca()

    # Main plot

    da = ds_le_percountry_GMT['mmm'].loc[{
        'country':country,
        'birth_year': np.arange(1960, 2021),
    }]

    max_val = da.max(skipna=True).item()

    if np.isnan(max_val):

        print("Plot not performed because NaN values are present in Lifetime Exposure for {}".format(country))
    
    else:

        if flags['extr']=='heatwavedarea':
        
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
    
            # Compute the maximum value of the data
            max_val = da.max(skipna=True).item()

            # Set a target and maximum number of contour levels
            target_n_levels = 30
            max_n_levels = 50

            # Compute raw step
            step = max_val / target_n_levels if target_n_levels != 0 else 0.1

            # Add a minimum threshold for the step to avoid too many levels
            min_step = 0.1
            rounded_step = max(np.round(step, 2), min_step)

            # Compute levels, but clip number of levels to avoid overflow
            n_levels = int(np.ceil(max_val / rounded_step)) + 1

            # Safety check to avoid exceeding the max allowed size
            if n_levels > max_n_levels:
                rounded_step = max_val / max_n_levels
                n_levels = max_n_levels

            # Create levels array
            levels_hw = np.arange(0, rounded_step * n_levels, rounded_step)

            # Safety check before plotting
            if len(np.unique(levels_hw)) < 2:

                print("Plot not performed: not enough contour levels for {}".format(country))
                return

            else:

                # Target and maximum number of ticks for the colorbar (for readability)
                target_n_ticks = 6
                max_n_ticks = 10

                # Compute a raw tick step
                tick_step = max_val / target_n_ticks if target_n_ticks != 0 else 0.5

                # Add a minimum threshold for the tick step
                min_tick_step = 0.1
                rounded_tick_step = max(np.round(tick_step, 2), min_tick_step)

                # Compute number of ticks, enforce a max limit
                n_ticks = int(np.ceil(max_val / rounded_tick_step)) + 1

                if n_ticks > max_n_ticks:
                    rounded_tick_step = max_val / max_n_ticks
                    n_ticks = max_n_ticks

                # Generate tick positions safely
                ticks_z = np.arange(0, rounded_tick_step * n_ticks, rounded_tick_step)

                ticks_z = np.round(ticks_z, 1)

                # Use in colorbar kwargs
                cbar_kwargs = {'ticks': ticks_z}

                p =da.plot.contourf(
                    x='birth_year',
                    y='GMT',
                    ax=ax,
                    add_labels=False,
                    levels=levels_hw,
                    cmap='Reds',
                    cbar_kwargs={'ticks': ticks_z}
                )

            ax.set_yticks(
                    ticks=GMT_indices_ticks,
                    labels=gmts2100,
                    color='gray',
            )
            
            ax.set_xticks(
                ticks=np.arange(1960,2025,10),
                color='gray',
            )  

            ax.tick_params(axis='both', labelsize=12)  

            ax.set_title(
                "{} : Lifetime Exposure in \n {}".format(flags['extr'], country),
                loc='center',
                fontweight='bold',
                color='gray',
                fontsize=16,
            )
            ax.spines['right'].set_color('gray')
            ax.spines['top'].set_color('gray')
            ax.spines['left'].set_color('gray')
            ax.spines['bottom'].set_color('gray')  
                
            ax.set_xlabel('Birth year',fontsize=14,color='gray')
            ax.set_ylabel('GMT warming by 2100 [°C]',fontsize=14,color='gray',rotation='vertical',)

            plt.tight_layout()
            
            plt.savefig(scripts_dir+'/figures/source2suffering/development/{}/fig8_BE_{}_lifetime_exposure_{}.png'.format(flags['extr'],flags['extr'],country))

            print("Plot performed for {}".format(country))

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

    # Plot lines with markers (instead of only points)
    plt.plot(birth_cohort_int, s2s_valc_nr_children_facing_extra_hazard, 
            '-o', color="tab:orange", label="Source2Suffering", markersize=10,lw=2)
    plt.plot(birth_cohort_int, wt_valc_nr_children_facing_extra_hazard, 
            '-o', color="tab:blue", label="W.Thiery", markersize=10,lw=2)

    # Set the figure title
    plt.title("Number of children facing at\nleast one extra {} due to NeptunDeep".format(hazards),
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
    plt.savefig(scripts_dir + '/figures/assessment/NeptunDeep/nr_children_facing_extra_{}_NeptunDeep_gmt_{}.png'.format(hazards, flags['gmt']))

    

