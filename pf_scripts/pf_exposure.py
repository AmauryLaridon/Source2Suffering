# --------------------------------------------------------------------------- #
# Functions to compute lifetime exposure to climate extreme                   #
# --------------------------------------------------------------------------- #

#%%---------------------------------------------------------------#
# Libraries                                                       #
# ----------------------------------------------------------------#

import sys
from operator import index
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import pickle as pk
from scipy import interpolate
from scipy import stats as sts
import regionmask as rm
import glob
import time
import matplotlib.pyplot as plt
from copy import deepcopy as cp
#%%---------------------------------------------------------------#
# Linear regression                                               #
#-----------------------------------------------------------------#

def lreg(x, y):
    # Wrapper around scipy linregress to use in apply_ufunc
    slope, intercept, r_value, p_value, std_err = sts.linregress(x, y)
    return np.array([slope, p_value, r_value])

#%%---------------------------------------------------------------#
# Apply vectorized linear regression                              #
#-----------------------------------------------------------------#

def vectorize_lreg(da_y, da_x=None):
    if da_x is not None:
        pass
    else:
        da_list = []
        for t in da_y.time.values:
            da_list.append(xr.where(da_y.sel(time=t).notnull(), t, da_y.sel(time=t)))
        da_x = xr.concat(da_list, dim='time')

    stats = xr.apply_ufunc(
        lreg,
        da_x,
        da_y,
        input_core_dims=[['time'], ['time']],
        output_core_dims=[["parameter"]],
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={"output_sizes": {"parameter": 3}},  
        output_dtypes=['float64']
    )

#%%---------------------------------------------------------------#
# Bootstrapping function                                          #
#-----------------------------------------------------------------#

def resample(
    da, 
    resample_dim,
    life_extent,
):
    """Resample with replacement in dimension ``resample_dim``. https://climpred.readthedocs.io/en/stable/_modules/climpred/bootstrap.html

    Args:
        initialized (xr.Dataset): input xr.Dataset to be resampled.
        resample_dim (str): dimension to resample along.
        life_extent (int): number of years per lifetime
        
    Returns:
        xr.Dataset: resampled along ``resample_dim``.

    """
    to_be_resampled = da[resample_dim].values
    smp = np.random.choice(to_be_resampled, life_extent)
    smp_da = da.sel({resample_dim: smp})
    smp_da[resample_dim] = np.arange(1960,1960+life_extent)
    return smp_da

#%%---------------------------------------------------------------#
# Improved function to compute extreme event exposure across a    #
# person's lifetime. Translation of W.Thiery's lifetime_exposure()#
# Matlab function done by L.Grant.                                #
#-----------------------------------------------------------------#

def calc_life_exposure(
    df_exposure,
    df_life_expectancy,
    col,
):

    # initialise birth years 
    exposure_birthyears_percountry = np.empty(len(df_life_expectancy))

    for i, birth_year in enumerate(df_life_expectancy.index):

        life_expectancy = df_life_expectancy.loc[birth_year,col] 

        # define death year based on life expectancy
        death_year = birth_year + np.floor(life_expectancy)

        # integrate exposure over full years lived
        exposure_birthyears_percountry[i] = df_exposure.loc[birth_year:death_year,col].sum()

        # add exposure during last (partial) year
        exposure_birthyears_percountry[i] = exposure_birthyears_percountry[i] + \
            df_exposure.loc[death_year+1,col].sum() * \
                (life_expectancy - np.floor(life_expectancy))

    # a series for each column to somehow group into a dataframe
    exposure_birthyears_percountry = pd.Series(
        exposure_birthyears_percountry,
        index=df_life_expectancy.index,
        name=col,
    )

    return exposure_birthyears_percountry

#%%---------------------------------------------------------------#
# Calculated weighted fieldmean per country mask                  #
#-----------------------------------------------------------------#

def calc_weighted_fldmean(
    da, 
    weights, 
    countries_mask,
    ind_country, 
    flag_region,
):

    # one country provided, easy masking
    # only keeps the AFA data for the country under study, for the others a NaN value is attributedppui
    if not flag_region : 
        da_masked = da.where(countries_mask == ind_country)

    # if more countries are provided, combine the different masks 
    else: 
        
        if len(ind_country) > 1:
            
            mask = xr.DataArray(
                np.in1d(countries_mask,ind_country).reshape(countries_mask.shape),
                dims=countries_mask.dims,
                coords=countries_mask.coords,
            )
            da_masked = da.where(mask)
    
    # weight the AFA of the country under study by the size of its population over time at the grid cell
    da_weighted_fldmean = da_masked.weighted(weights).mean(dim=("lat", "lon"))

    return da_weighted_fldmean

#%%---------------------------------------------------------------#
# Get member countries per region                                 #
#-----------------------------------------------------------------#

def get_countries_of_region(
    region, 
    df_countries,
): 

    # Get list of member countries from region
    member_countries = df_countries.loc[df_countries['region']==region]['name'].values

    # not region but income group
    if len(member_countries) == 0: 
        member_countries = df_countries.loc[df_countries['incomegroup']==region]['name'].values

    # get all countries for the world
    if region == 'World':
        member_countries = df_countries['name'].values

    return member_countries    

#%%---------------------------------------------------------------#
# Function to compute multi-model mean across ISIMIP simulations  #
# based on mf_exposure_mmm.m (see Thiery et al.(2021))            #
#-----------------------------------------------------------------#

def calc_lifetime_exposure_mmm_xr(
    ds_le_perrun,
    flags,
):
    # remove from the output a warning that appears often  
    import warnings
    warnings.filterwarnings("ignore", message="All-NaN slice encountered")

    if 'country' in ds_le_perrun.dims:

        # ------------------------ OverShoot (OS) trajectory ------------------------- #
        
        mmm_OS = ds_le_perrun['le_percountry_perrun_OS'].mean(dim='run', skipna=True)
        std_OS = ds_le_perrun['le_percountry_perrun_OS'].std(dim='run', skipna=True)
        lqntl_OS = ds_le_perrun['le_percountry_perrun_OS'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_OS = ds_le_perrun['le_percountry_perrun_OS'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_OS = ds_le_perrun['le_percountry_perrun_OS'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_OS = mmm_OS / mmm_OS.sel(birth_year=1960)
        lqntl_EMF_OS = lqntl_OS / mmm_OS.sel(birth_year=1960)
        uqntl_EMF_OS = uqntl_OS / mmm_OS.sel(birth_year=1960)
        median_EMF_OS = median_OS / mmm_OS.sel(birth_year=1960)

        # --------------------- no-OverShoot (noOS) trajectory ---------------------- #

        mmm_noOS = ds_le_perrun['le_percountry_perrun_noOS'].mean(dim='run', skipna=True)
        std_noOS = ds_le_perrun['le_percountry_perrun_noOS'].std(dim='run', skipna=True)
        lqntl_noOS = ds_le_perrun['le_percountry_perrun_noOS'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_noOS = ds_le_perrun['le_percountry_perrun_noOS'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_noOS = ds_le_perrun['le_percountry_perrun_noOS'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_noOS = mmm_noOS / mmm_noOS.sel(birth_year=1960)
        lqntl_EMF_noOS = lqntl_noOS / mmm_noOS.sel(birth_year=1960)
        uqntl_EMF_noOS = uqntl_noOS / mmm_noOS.sel(birth_year=1960)
        median_EMF_noOS = median_noOS / mmm_noOS.sel(birth_year=1960)

        # ------------------------ STS-ModAct trajectory ------------------------- #
        
        mmm_STS_ModAct = ds_le_perrun['le_percountry_perrun_STS_ModAct'].mean(dim='run', skipna=True)
        std_STS_ModAct = ds_le_perrun['le_percountry_perrun_STS_ModAct'].std(dim='run', skipna=True)
        lqntl_STS_ModAct = ds_le_perrun['le_percountry_perrun_STS_ModAct'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_STS_ModAct = ds_le_perrun['le_percountry_perrun_STS_ModAct'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_STS_ModAct = ds_le_perrun['le_percountry_perrun_STS_ModAct'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_STS_ModAct = mmm_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)
        lqntl_EMF_STS_ModAct = lqntl_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)
        uqntl_EMF_STS_ModAct = uqntl_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)
        median_EMF_STS_ModAct = median_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)

        # ------------------------ STS-Ren trajectory ------------------------- #
        
        mmm_STS_Ren = ds_le_perrun['le_percountry_perrun_STS_Ren'].mean(dim='run', skipna=True)
        std_STS_Ren = ds_le_perrun['le_percountry_perrun_STS_Ren'].std(dim='run', skipna=True)
        lqntl_STS_Ren = ds_le_perrun['le_percountry_perrun_STS_Ren'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_STS_Ren = ds_le_perrun['le_percountry_perrun_STS_Ren'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_STS_Ren = ds_le_perrun['le_percountry_perrun_STS_Ren'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_STS_Ren = mmm_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)
        lqntl_EMF_STS_Ren = lqntl_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)
        uqntl_EMF_STS_Ren = uqntl_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)
        median_EMF_STS_Ren = median_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)

        # --------------------- Stylized trajecotries of the BE ---------------------- #

        mmm_BE = ds_le_perrun['le_percountry_perrun_BE'].mean(dim='run', skipna=True)
        std_BE = ds_le_perrun['le_percountry_perrun_BE'].std(dim='run', skipna=True)
        lqntl_BE = ds_le_perrun['le_percountry_perrun_BE'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_BE = ds_le_perrun['le_percountry_perrun_BE'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_BE = ds_le_perrun['le_percountry_perrun_BE'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_BE = mmm_BE / mmm_BE.sel(birth_year=1960)
        lqntl_EMF_BE = lqntl_BE / mmm_BE.sel(birth_year=1960)
        uqntl_EMF_BE = uqntl_BE / mmm_BE.sel(birth_year=1960)
        median_EMF_BE = median_BE / mmm_BE.sel(birth_year=1960)

    if 'region' in ds_le_perrun.dims:
    
        # ------------------------ OverShoot (OS) trajectory ------------------------- #
        
        mmm_OS = ds_le_perrun['le_perregion_perrun_OS'].mean(dim='run', skipna=True)
        std_OS = ds_le_perrun['le_perregion_perrun_OS'].std(dim='run', skipna=True)
        lqntl_OS = ds_le_perrun['le_perregion_perrun_OS'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_OS = ds_le_perrun['le_perregion_perrun_OS'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        median_OS = ds_le_perrun['le_perregion_perrun_OS'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_OS = mmm_OS / mmm_OS.sel(birth_year=1960)
        lqntl_EMF_OS = lqntl_OS / mmm_OS.sel(birth_year=1960)
        uqntl_EMF_OS = uqntl_OS / mmm_OS.sel(birth_year=1960)
        median_EMF_OS = median_OS / mmm_OS.sel(birth_year=1960)

        # --------------------- no-OverShoot (noOS) trajectory ---------------------- #

        mmm_noOS = ds_le_perrun['le_perregion_perrun_noOS'].mean(dim='run', skipna=True)
        std_noOS = ds_le_perrun['le_perregion_perrun_noOS'].std(dim='run', skipna=True)
        lqntl_noOS = ds_le_perrun['le_perregion_perrun_noOS'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_noOS = ds_le_perrun['le_perregion_perrun_noOS'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_noOS = ds_le_perrun['le_perregion_perrun_noOS'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_noOS = mmm_noOS / mmm_noOS.sel(birth_year=1960)
        lqntl_EMF_noOS = lqntl_noOS / mmm_noOS.sel(birth_year=1960)
        uqntl_EMF_noOS = uqntl_noOS / mmm_noOS.sel(birth_year=1960)
        median_EMF_noOS = median_noOS / mmm_noOS.sel(birth_year=1960)

        # ------------------------ STS-ModAct trajectory ------------------------- #
        
        mmm_STS_ModAct = ds_le_perrun['le_perregion_perrun_STS_ModAct'].mean(dim='run', skipna=True)
        std_STS_ModAct = ds_le_perrun['le_perregion_perrun_STS_ModAct'].std(dim='run', skipna=True)
        lqntl_STS_ModAct = ds_le_perrun['le_perregion_perrun_STS_ModAct'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_STS_ModAct = ds_le_perrun['le_perregion_perrun_STS_ModAct'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        median_STS_ModAct = ds_le_perrun['le_perregion_perrun_STS_ModAct'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_STS_ModAct = mmm_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)
        lqntl_EMF_STS_ModAct = lqntl_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)
        uqntl_EMF_STS_ModAct = uqntl_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)
        median_EMF_STS_ModAct = median_STS_ModAct / mmm_STS_ModAct.sel(birth_year=1960)

        # ------------------------ STS-Ren trajectory ------------------------- #
        
        mmm_STS_Ren = ds_le_perrun['le_perregion_perrun_STS_Ren'].mean(dim='run', skipna=True)
        std_STS_Ren = ds_le_perrun['le_perregion_perrun_STS_Ren'].std(dim='run', skipna=True)
        lqntl_STS_Ren = ds_le_perrun['le_perregion_perrun_STS_Ren'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_STS_Ren = ds_le_perrun['le_perregion_perrun_STS_Ren'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        median_STS_Ren = ds_le_perrun['le_perregion_perrun_STS_Ren'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_STS_Ren = mmm_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)
        lqntl_EMF_STS_Ren = lqntl_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)
        uqntl_EMF_STS_Ren = uqntl_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)
        median_EMF_STS_Ren = median_STS_Ren / mmm_STS_Ren.sel(birth_year=1960)

        # --------------------- Stylized trajecotries of the BE ---------------------- #

        mmm_BE = ds_le_perrun['le_perregion_perrun_BE'].mean(dim='run', skipna=True)
        std_BE = ds_le_perrun['le_perregion_perrun_BE'].std(dim='run', skipna=True)
        lqntl_BE = ds_le_perrun['le_perregion_perrun_BE'].quantile(
            q=0.25,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        uqntl_BE = ds_le_perrun['le_perregion_perrun_BE'].quantile(
            q=0.75,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )

        median_BE = ds_le_perrun['le_perregion_perrun_BE'].quantile(
            q=0.5,
            dim='run',
            method='inverted_cdf',
            skipna=True
        )
        
        # get EMF of stats (divide by 60 yr old / 1960 cohort)
        mmm_EMF_BE = mmm_BE / mmm_BE.sel(birth_year=1960)
        lqntl_EMF_BE = lqntl_BE / mmm_BE.sel(birth_year=1960)
        uqntl_EMF_BE = uqntl_BE / mmm_BE.sel(birth_year=1960)
        median_EMF_BE = median_BE / mmm_BE.sel(birth_year=1960)

    # --------------------- Integration in the DataSet ---------------------------- #

    ds_le_perrun['mmm_OS'] = mmm_OS
    ds_le_perrun['std_OS'] = std_OS
    ds_le_perrun['lqntl_OS'] = lqntl_OS
    ds_le_perrun['uqntl_OS'] = uqntl_OS
    ds_le_perrun['median_OS'] = median_OS
    ds_le_perrun['mmm_EMF_OS'] = mmm_EMF_OS
    ds_le_perrun['lqntl_EMF_OS'] = lqntl_EMF_OS
    ds_le_perrun['uqntl_EMF_OS'] = uqntl_EMF_OS
    ds_le_perrun['median_EMF_OS'] = median_EMF_OS 

    ds_le_perrun['mmm_noOS'] = mmm_noOS
    ds_le_perrun['std_noOS'] = std_noOS
    ds_le_perrun['lqntl_noOS'] = lqntl_noOS
    ds_le_perrun['uqntl_noOS'] = uqntl_noOS
    ds_le_perrun['median_noOS'] = median_noOS
    ds_le_perrun['mmm_EMF_noOS'] = mmm_EMF_noOS
    ds_le_perrun['lqntl_EMF_noOS'] = lqntl_EMF_noOS
    ds_le_perrun['uqntl_EMF_noOS'] = uqntl_EMF_noOS
    ds_le_perrun['median_EMF_noOS'] = median_EMF_noOS  

    ds_le_perrun['mmm_STS_ModAct'] = mmm_STS_ModAct
    ds_le_perrun['std_STS_ModAct'] = std_STS_ModAct
    ds_le_perrun['lqntl_STS_ModAct'] = lqntl_STS_ModAct
    ds_le_perrun['uqntl_STS_ModAct'] = uqntl_STS_ModAct
    ds_le_perrun['median_STS_ModAct'] = median_STS_ModAct
    ds_le_perrun['mmm_EMF_STS_ModAct'] = mmm_EMF_STS_ModAct
    ds_le_perrun['lqntl_EMF_STS_ModAct'] = lqntl_EMF_STS_ModAct
    ds_le_perrun['uqntl_EMF_STS_ModAct'] = uqntl_EMF_STS_ModAct
    ds_le_perrun['median_EMF_STS_ModAct'] = median_EMF_STS_ModAct   

    ds_le_perrun['mmm_STS_Ren'] = mmm_STS_Ren
    ds_le_perrun['std_STS_Ren'] = std_STS_Ren
    ds_le_perrun['lqntl_STS_Ren'] = lqntl_STS_Ren
    ds_le_perrun['uqntl_STS_Ren'] = uqntl_STS_Ren
    ds_le_perrun['median_STS_Ren'] = median_STS_Ren
    ds_le_perrun['mmm_EMF_STS_Ren'] = mmm_EMF_STS_Ren
    ds_le_perrun['lqntl_EMF_STS_Ren'] = lqntl_EMF_STS_Ren
    ds_le_perrun['uqntl_EMF_STS_Ren'] = uqntl_EMF_STS_Ren
    ds_le_perrun['median_EMF_STS_Ren'] = median_EMF_STS_Ren

    ds_le_perrun['mmm_BE'] = mmm_BE
    ds_le_perrun['std_BE'] = std_BE
    ds_le_perrun['lqntl_BE'] = lqntl_BE
    ds_le_perrun['uqntl_BE'] = uqntl_BE
    ds_le_perrun['median_BE'] = median_BE
    ds_le_perrun['mmm_EMF_BE'] = mmm_EMF_BE
    ds_le_perrun['lqntl_EMF_BE'] = lqntl_EMF_BE
    ds_le_perrun['uqntl_EMF_BE'] = uqntl_EMF_BE
    ds_le_perrun['median_EMF_BE'] = median_EMF_BE 
    

    if 'region' in ds_le_perrun.dims:

        # dump pickle of lifetime exposure per region
        with open(data_dir+'{}/{}/ds_le_perregion_gmt_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
            pk.dump(ds_le_perrun,f)
    
    if 'country' in ds_le_perrun.dims:

        # dump pickle of lifetime exposure per country
        with open(data_dir+'{}/{}/ds_le_percountry_gmt_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
            pk.dump(ds_le_perrun,f)

    print("\nEnd of the computation\n")

    return ds_le_perrun
    
#%%---------------------------------------------------------------#
# Function to compute multi-model mean across ISIMIP simulations  #
# based on mf_exposure_mmm.m (see Thiery et al.(2021))            #
#-----------------------------------------------------------------#

def calc_exposure_mmm_pic_xr(
    d_exposure_pic,
    dim_1_name,
    var_tag,
):        
        
    # concat pic data array from dict of separate arrays
    da_exposure_pic = xr.concat(
        [v for v in d_exposure_pic.values()],
        dim='runs',    
    ).assign_coords({'runs':list(d_exposure_pic.keys())})

    # runs and lifetimes (lifetimes from boostrapping) redundant, so compile together
    da_exposure_pic = da_exposure_pic.stack(
        pic_lifetimes=['runs','lifetimes'],
    )

    # pic exposure stats for EMF
    da_exposure_pic_mmm = da_exposure_pic.mean(dim='pic_lifetimes')
    da_exposure_pic_std = da_exposure_pic.std(dim='pic_lifetimes')
    da_exposure_pic_lqntl = da_exposure_pic.quantile(
        q=0.25,
        dim='pic_lifetimes',
        method='inverted_cdf',
    )
    da_exposure_pic_uqntl = da_exposure_pic.quantile(
        q=0.75,
        dim='pic_lifetimes',
        method='inverted_cdf',
    )    

    # pic quantile for birth cohort exposure emergence
    da_exposure_pic_ext = da_exposure_pic.quantile(
        q=0.9999,
        dim='pic_lifetimes',
        # method='inverted_cdf',
    )
    
    # assemble into dataset
    ds_exposure_pic_stats = xr.Dataset(
        data_vars={
            'mmm_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_mmm.data),
            'std_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_std.data),
            'lqntl_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_lqntl.data),
            'uqntl_{}'.format(var_tag): ([dim_1_name],da_exposure_pic_uqntl.data),
            'ext': ([dim_1_name],da_exposure_pic_ext.data),
        },
        coords={
            dim_1_name: (dim_1_name,da_exposure_pic[dim_1_name].data),
        }
    )

    return ds_exposure_pic_stats

#%%---------------------------------------------------------------#
# Convert Area Fraction Affected (AFA) to per-country and bassins #
# number of extremes affecting one individual across life span.   #
# Compute also the trend of these exposure                        #
#-----------------------------------------------------------------#

def calc_exposure_trends(
    d_isimip_meta,
    grid_area,
    gdf_country_borders,
    flags,
):

    # time step in year used for the assessment
    time_res = 1

    # arrays of lat/lon values
    lat = grid_area.lat.values
    lon = grid_area.lon.values
    
    # 3d mask for ar6 regions
    ar6_regs_3D = rm.defined_regions.ar6.land.mask_3D(lon,lat)
    
    # 3d mask for countries
    countries_3D = rm.mask_3D_geopandas(gdf_country_borders.reset_index(),lon,lat)
    
    # basin shapefiles
    gdf_basins = gpd.read_file(data_dir+'shapefiles/Major_Basins_of_the_World.shp')
    gdf_basins = gdf_basins.loc[:,['NAME','geometry']]
    
    # merge basins with multiple entries
    basins_grouped = []
    bc = {k:0 for k in gdf_basins['NAME']} # bc for basin counter
    for b_name in gdf_basins['NAME']:
        if len(gdf_basins.loc[gdf_basins['NAME']==b_name]) > 1:
            if bc[b_name]==0:
                gdf_basin = gdf_basins.loc[gdf_basins['NAME']==b_name]
                basins_grouped.append(gdf_basin.dissolve())
            bc[b_name]+=1
        else:
            basins_grouped.append(gdf_basins.loc[gdf_basins['NAME']==b_name])
    gdf_basins = pd.concat(basins_grouped).reset_index().loc[:,['NAME','geometry']]
    basins_3D = rm.mask_3D_geopandas(gdf_basins,lon,lat) # 3d mask for basins

    # dataset for exposure trends
    ds_le_trends_regions = xr.Dataset(
        data_vars={
            'exposure_trend_ar6': (
                ['run','GMT','region','year'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(GMT_labels),len(ar6_regs_3D.region.data),len(np.arange(year_start,year_ref+1,time_res))),
                    fill_value=np.nan,
                ),
            ),
            'mean_exposure_trend_ar6': (
                ['GMT','region','year'],
                np.full(
                    (len(GMT_labels),len(ar6_regs_3D.region.data),len(np.arange(year_start,year_ref+1,time_res))),
                    fill_value=np.nan,
                ),
            ),                                
            'exposure_trend_country': (
                ['run','GMT','country','year'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(GMT_labels),len(countries_3D.region.data),len(np.arange(year_start,year_ref+1,time_res))),
                    fill_value=np.nan,
                ),
            ),
            'mean_exposure_trend_country': (
                ['GMT','country','year'],
                np.full(
                    (len(GMT_labels),len(countries_3D.region.data),len(np.arange(year_start,year_ref+1,time_res))),
                    fill_value=np.nan,
                ),
            ),         
            'exposure_trend_basin': (
                ['run','GMT','basin','year'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(GMT_labels),len(basins_3D.region.data),len(np.arange(year_start,year_ref+1,time_res))),
                    fill_value=np.nan,
                ),
            ),
            'mean_exposure_trend_basin': (
                ['GMT','basin','year'],
                np.full(
                    (len(GMT_labels),len(basins_3D.region.data),len(np.arange(year_start,year_ref+1,time_res))),
                    fill_value=np.nan,
                ),
            ),                                                             
        },
        coords={
            'region': ('region', ar6_regs_3D.region.data),
            'country': ('country', countries_3D.region.data),
            'basin': ('basin', basins_3D.region.data),
            'run': ('run', list(d_isimip_meta.keys())),
            'GMT': ('GMT', GMT_labels),
            'year': ('year', np.arange(year_start,year_ref+1,time_res))
        }
    )
    
    # loop over simulations
    for i in list(d_isimip_meta.keys()): 

        print('ISIMIP Simulation {} of {}'.format(i,len(d_isimip_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA = pk.load(f)  
        
        # per GMT step, if max threshold criteria met, run gmt mapping and compute trends
        for step in GMT_labels:
                
            if d_isimip_meta[i]['GMT_strj_valid'][step]:

                print("Simulations used for re-mapping for GMT index = {}".format([step]))

                # reindexing original exposure array based on GMT-mapping indices
                da_AFA_step = da_AFA.reindex(
                    {'time':da_AFA['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
                ).assign_coords({'time':year_range}) 
                
                # get sums of exposed area per ar6, country & basin, convert m^2 to km^2
                da_AFA_ar6_weighted_sum = da_AFA_step.weighted(ar6_regs_3D*grid_area/10**6).sum(dim=('lat','lon'))
                da_AFA_country_weighted_sum = da_AFA_step.weighted(countries_3D*grid_area/10**6).sum(dim=('lat','lon'))
                da_AFA_basin_weighted_sum = da_AFA_step.weighted(basins_3D*grid_area/10**6).sum(dim=('lat','lon'))
        
                # regressions on separate 80 year periods of area sums
                for y in np.arange(year_start,year_ref+1,time_res):
                    
                    # ar6 regions
                    ds_le_trends_regions['exposure_trend_ar6'].loc[{
                        'run':i,
                        'GMT':step,
                        'region':ar6_regs_3D.region.data,
                        'year':y,
                    }] = vectorize_lreg(da_AFA_ar6_weighted_sum.loc[{'time':np.arange(y,y+81)}])
                    
                    # countries
                    ds_le_trends_regions['exposure_trend_country'].loc[{
                        'run':i,
                        'GMT':step,
                        'country':countries_3D.region.data,
                        'year':y,
                    }] = vectorize_lreg(da_AFA_country_weighted_sum.loc[{'time':np.arange(y,y+81)}])
                    
                    # basins
                    ds_le_trends_regions['exposure_trend_basin'].loc[{
                        'run':i,
                        'GMT':step,
                        'basin':basins_3D.region.data,
                        'year':y,
                    }] = vectorize_lreg(da_AFA_basin_weighted_sum.loc[{'time':np.arange(y,y+81)}])
    
    # take means of trends in exposed area
    ds_le_trends_regions['mean_exposure_trend_ar6'] = ds_le_trends_regions['exposure_trend_ar6'].mean(dim='run')
    ds_le_trends_regions['mean_exposure_trend_country'] = ds_le_trends_regions['exposure_trend_country'].mean(dim='run')
    ds_le_trends_regions['mean_exposure_trend_basin'] = ds_le_trends_regions['exposure_trend_basin'].mean(dim='run')

    # dump pickle of exposure trends
    with open(data_dir+'{}/{}/lifetime_exposure_trends_regions.pkl'.format(flags['version'],flags['extr']), 'wb') as f:
        pk.dump(ds_le_trends_regions,f)

    return ds_le_trends_regions
        
#%%----------------------------------------------------------------#
# Convert Area Fraction Affected (AFA) to per-country number of    #
# extremes affecting one individual across life span               #
# Original W.Thiery's lifetime exposure function translated by     # 
# L.Grant even thought not used for Grant et al.(2025) since the   #
# lifetime exposure is computed in the emergence analysis          #
# This function is used for the Source2Suffering and Laridon       #
# et al.(2025) analysis for backward compatibility with W.Thiery's #
# results. A.Laridon have add the exposure per region per run      # 
# computation following Thiery et al.(2021)                        #
#----------------------------------------------------------------- #

def calc_lifetime_exposure(
    d_isimip_meta, 
    df_countries, 
    countries_regions, 
    countries_mask, 
    da_population, 
    df_life_expectancy_5,
    ds_regions,
    d_cohort_weights_regions,
    flags,
):


    #---------------------------------------------------------------------#
    # Init                                                                #
    #---------------------------------------------------------------------#

    # nan dataset of lifetime exposure at the country scale
    ds_le_percountry_perrun_GMT = xr.Dataset(
        data_vars={
            'lifetime_exposure': (
                ['run','GMT','country','birth_year'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(GMT_labels),len(df_countries['name'].values),len(birth_years)),
                    fill_value=np.nan,
                ),
            )
        },
        coords={
            'country': ('country', df_countries['name'].values),
            'birth_year': ('birth_year', birth_years),
            'run': ('run', np.arange(1,len(list(d_isimip_meta.keys()))+1)),
            'GMT': ('GMT', GMT_labels)
        }
    )

    nregions = len(ds_regions['name']) # number of regions that will be used for the init of the dimension

    # nan dataset of lifetime exposure at the country scale
    ds_le_perregion_perrun_GMT = xr.Dataset(
        data_vars={
            'lifetime_exposure': (
                ['run','GMT','region','birth_year'],
                np.full(
                    (len(list(d_isimip_meta.keys())),len(GMT_labels),nregions,len(birth_years)),
                    fill_value=np.nan,
                ),
            )
        },
        coords={
            'region': ('region', np.arange(0,nregions,1)),
            'birth_year': ('birth_year', birth_years),
            'run': ('run', np.arange(1,len(list(d_isimip_meta.keys()))+1)),
            'GMT': ('GMT', GMT_labels)
        }
    )
    
    #---------------------------------------------------------------------#
    # Loop over ISIMIP simulations                                        #
    #---------------------------------------------------------------------#
    for i in list(d_isimip_meta.keys()): 

        print('ISIMIP Simulation {} of {}'.format(i,len(d_isimip_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA = pk.load(f)

        #---------------------------------------------------------------------#
        # Per country                                                         #
        #---------------------------------------------------------------------#

        # initialise dicts
        d_exposure_peryear_percountry = {}

        # get spatial average
        for j, country in enumerate(df_countries['name']):

            print('Processing country '+str(j+1)+' of '+str(len(df_countries)), end='\r')
            
            # calculate mean per country weighted by population
            ind_country = countries_regions.map_keys(country)

            # historical + RCP simulations
            d_exposure_peryear_percountry[country] = calc_weighted_fldmean( 
                da_AFA,
                da_population, 
                countries_mask, 
                ind_country, 
                flag_region=False,
            )
                        
        #---------------------------------------------------------------------#
        # Convert dict to dataframe for vectorizing and integrate exposures   # 
        #---------------------------------------------------------------------#
       
        frame = {k:v.values for k,v in d_exposure_peryear_percountry.items()}
        df_exposure = pd.DataFrame(frame,index=year_range)           

        #---------------------------------------------------------------------#
        # Computations for Burning Embers diagram by mapping the GMTs         #
        #---------------------------------------------------------------------#   
             
        for step in GMT_labels:

            # if max threshold criteria met, run gmt mapping
            if d_isimip_meta[i]['GMT_strj_valid'][step]:

                print("Simulations used for re-mapping for GMT index = {}".format([step]))
            
                # reindexing original exposure array based on GMT-mapping indices
                d_exposure_perrun_step = df_exposure.apply(
                    lambda col: calc_life_exposure(
                        df_exposure.reindex(df_exposure.index[d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]).set_index(df_exposure.index),
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
        
                # convert dataframe to data array of lifetime exposure (le) per country and birth year
                ds_le_percountry_perrun_GMT['lifetime_exposure'].loc[{
                    'run':i,
                    'GMT':step,
                }] = d_exposure_perrun_step.values.transpose() 
    
                #---------------------------------------------------------------------#
                # Per region                                                          #
                #---------------------------------------------------------------------#

                for region_ind, region in enumerate(ds_regions.region.values):
                
                    region_name = ds_regions['name'].sel(region=region_ind).item()
                    member_countries = ds_regions['member_countries'].sel(region=region_ind).values.tolist()

                    weights = d_cohort_weights_regions[region_name]

                    le_percountry = ds_le_percountry_perrun_GMT['lifetime_exposure'].loc[{
                        'run': i,
                        'GMT': step,
                        'country': member_countries,
                        'birth_year': slice(None)
                    }]

                    # print('region name', region_name)

                    # print('le_percountry')
                    # print(le_percountry)
                    # print(type(le_percountry))
                    # print(np.shape(le_percountry))

                    # print('weights')
                    # print(weights)
                    # print(type(weights))
                    # print(np.shape(weights))

                    weighted_avg = np.nansum(le_percountry.values * weights.values.T, axis=0) / np.nansum(weights.values.T, axis=0)

                    # print('weighted_avg')
                    # print(weighted_avg)
                    # print(type(weighted_avg))
                    # print(np.shape(weighted_avg))

                    ds_le_perregion_perrun_GMT['lifetime_exposure'].loc[{
                        'run': i,
                        'GMT': step,
                        'region': region_ind,
                        'birth_year': slice(None)
                    }] = weighted_avg


    # dump pickle of lifetime exposure per country
    with open(data_dir+'{}/{}/ds_le_percountry_perrun_{}_GMT.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
        pk.dump(ds_le_percountry_perrun_GMT,f)
    
    # dump pickle of lifetime exposure per region
    with open(data_dir+'{}/{}/ds_le_perregion_perrun_{}_GMT.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
        pk.dump(ds_le_perregion_perrun_GMT,f)

    return ds_le_percountry_perrun_GMT, ds_le_perregion_perrun_GMT

def calc_lifetime_exposure_v2(
    d_isimip_meta, 
    df_countries, 
    countries_regions, 
    countries_mask, 
    da_population, 
    df_life_expectancy_5,
    ds_regions,
    d_cohort_weights_regions,
    flags,
):


    #---------------------------------------------------------------------#
    # Init                                                                #
    #---------------------------------------------------------------------#

    ### Init of the ds_le_percountry_perrun DataSet ### 

    ds_le_percountry_perrun = xr.Dataset(

    data_vars={

        'le_percountry_perrun_OS': (
            ['run', 'country', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    len(df_countries['name'].values),
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_percountry_perrun_noOS': (
            ['run', 'country', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    len(df_countries['name'].values),
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_percountry_perrun_STS_ModAct': (
            ['run', 'country', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    len(df_countries['name'].values),
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_percountry_perrun_STS_Ren': (
            ['run', 'country', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    len(df_countries['name'].values),
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_percountry_perrun_BE': (
            ['run', 'country', 'birth_year', 'GMT'],
            np.full(
                (
                    len(d_isimip_meta),
                    len(df_countries['name'].values),
                    len(birth_years),
                    len(GMT_labels)
                ),
                fill_value=np.nan
            )
        )
    },

    coords={
        'run': ('run', np.arange(1, len(d_isimip_meta) + 1)),
        'country': ('country', df_countries['name'].values),
        'birth_year': ('birth_year', birth_years),
        'GMT': ('GMT', GMT_labels)
    }
    )

    nregions = len(ds_regions['name']) # number of regions that will be used for the init of the dimension

    ### Init of the ds_le_perregion_perrun DataSet ###

    ds_le_perregion_perrun = xr.Dataset(

    data_vars={

        'le_perregion_perrun_OS': (
            ['run', 'region', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    nregions,
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_perregion_perrun_noOS': (
            ['run', 'region', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    nregions,
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_perregion_perrun_STS_ModAct': (
            ['run', 'region', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    nregions,
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_perregion_perrun_STS_Ren': (
            ['run', 'region', 'birth_year'],
            np.full(
                (
                    len(d_isimip_meta),
                    nregions,
                    len(birth_years)
                ),
                fill_value=np.nan
            )
        ),

        'le_perregion_perrun_BE': (
            ['run', 'region', 'birth_year', 'GMT'],
            np.full(
                (
                    len(d_isimip_meta),
                    nregions,
                    len(birth_years),
                    len(GMT_labels)
                ),
                fill_value=np.nan
            )
        )
    },

    coords={
        'run': ('run', np.arange(1, len(d_isimip_meta) + 1)),
        'region': ('region', np.arange(0, nregions)),
        'birth_year': ('birth_year', birth_years),
        'GMT': ('GMT', GMT_labels)
    }
    )
    
    #---------------------------------------------------------------------#
    # Loop over ISIMIP simulations                                        #
    #---------------------------------------------------------------------#
    for i in list(d_isimip_meta.keys()): 

        print('\n---------- ISIMIP Simulation {} of {} ----------'.format(i,len(d_isimip_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA = pk.load(f)

        #---------------------------------------------------------------------#
        # Computation of the weighted field mean of AFA for each ISIMIP       #
        # simulations for each country                                        #
        #---------------------------------------------------------------------#

        # initialise dicts
        d_exposure_percountry = {}

        # get spatial average
        for j, country in enumerate(df_countries['name']):

            print('Computing the Spatial Average of the Exposure for country '+str(j+1)+' of '+str(len(df_countries)), end='\r')
            
            # calculate mean per country weighted by population
            ind_country = countries_regions.map_keys(country)

            # historical + RCP simulations
            d_exposure_percountry[country] = calc_weighted_fldmean( 
                da_AFA,
                da_population, 
                countries_mask, 
                ind_country, 
                flag_region=False,
            )
                        
        # Convert dict to dataframe for vectorizing and integrate exposures   
       
        frame = {k:v.values for k,v in d_exposure_percountry.items()}
        df_exposure_percountry = pd.DataFrame(frame,index=year_range)         


        # ------------------------------------------------------------------- #
        # Computation of the Lifetime Exposure per country for the pre-design #
        # trajectories by mapping the GMTs of the ISMIP simulations to the    #
        # pre-design trajectories if the ISMIP simulations use is valid for   #
        # remapping                                                           #
        # ------------------------------------------------------------------- #

        #----- Computation of the Lifetime Exposure for the OverShoot (OS) Trajectory -----#

        print("                                                                          ")
        print('\nComputation of the Lifetime Exposure for the OverShoot (OS) trajectory')

        # if max threshold criteria met, run gmt mapping
        if d_isimip_meta[i]['GMT_OS_valid']:

            print("ISMIP Simulation {} use for re-mapping GMT index for the OS trajectory".format(i))

            # reindexing original exposure array based on GMT-mapping indices
            d_le_percountry_perrun_OS = df_exposure_percountry.apply(
                lambda col: calc_life_exposure(
                    df_exposure_percountry.reindex(df_exposure_percountry.index[d_isimip_meta[i]['ind_RCP2GMT_OS'][:]]).set_index(df_exposure_percountry.index),
                    df_life_expectancy_5,
                    col.name,
                ),
                axis=0,
            )

            # convert dataframe to data array of lifetime exposure (le) per country and birth year
            ds_le_percountry_perrun['le_percountry_perrun_OS'].loc[{
                'run':i,
            }] = d_le_percountry_perrun_OS.values.transpose() 

            #---------------------------------------------------------------------#
            # Per region                                                          #
            #---------------------------------------------------------------------#

            for region_ind, region in enumerate(ds_regions.region.values):
            
                region_name = ds_regions['name'].sel(region=region_ind).item()
                member_countries = ds_regions['member_countries'].sel(region=region_ind).values.tolist()

                weights = d_cohort_weights_regions[region_name]

                le_percountry_perrun = ds_le_percountry_perrun['le_percountry_perrun_OS'].loc[{
                    'run': i,
                    'country': member_countries,
                    'birth_year': slice(None)
                }]

                weighted_avg = np.nansum(le_percountry_perrun.values * weights.values.T, axis=0) / np.nansum(weights.values.T, axis=0)

                ds_le_perregion_perrun['le_perregion_perrun_OS'].loc[{
                    'run': i,
                    'region': region_ind,
                    'birth_year': slice(None)
                }] = weighted_avg
        
        else:

            print("ISMIP Simulation {} can not be used for re-mapping GMT for the OS trajectory".format(i))
        
        #----- Computation of the Lifetime Exposure for the no-OverShoot (noOS) Trajectory -----#

        print('\nComputation of the Lifetime Exposure for the no-OverShoot (noOS) trajectory')

        # if max threshold criteria met, run gmt mapping
        if d_isimip_meta[i]['GMT_noOS_valid']:

            print("ISMIP Simulation {} use for re-mapping GMT index for the noOS trajectory".format(i))

            # reindexing original exposure array based on GMT-mapping indices
            d_le_percountry_perrun_noOS = df_exposure_percountry.apply(
                lambda col: calc_life_exposure(
                    df_exposure_percountry.reindex(df_exposure_percountry.index[d_isimip_meta[i]['ind_RCP2GMT_noOS'][:]]).set_index(df_exposure_percountry.index),
                    df_life_expectancy_5,
                    col.name,
                ),
                axis=0,
            )

            # convert dataframe to data array of lifetime exposure (le) per country and birth year
            ds_le_percountry_perrun['le_percountry_perrun_noOS'].loc[{
                'run':i,
            }] = d_le_percountry_perrun_noOS.values.transpose() 

            #---------------------------------------------------------------------#
            # Per region                                                          #
            #---------------------------------------------------------------------#

            for region_ind, region in enumerate(ds_regions.region.values):
            
                region_name = ds_regions['name'].sel(region=region_ind).item()
                member_countries = ds_regions['member_countries'].sel(region=region_ind).values.tolist()

                weights = d_cohort_weights_regions[region_name]

                le_percountry_perrun = ds_le_percountry_perrun['le_percountry_perrun_noOS'].loc[{
                    'run': i,
                    'country': member_countries,
                    'birth_year': slice(None)
                }]

                weighted_avg = np.nansum(le_percountry_perrun.values * weights.values.T, axis=0) / np.nansum(weights.values.T, axis=0)

                ds_le_perregion_perrun['le_perregion_perrun_noOS'].loc[{
                    'run': i,
                    'region': region_ind,
                    'birth_year': slice(None)
                }] = weighted_avg
        
        else:

            print("ISMIP Simulation {} can not be used for re-mapping GMT for the noOS trajectory".format(i))

        #----- Computation of the Lifetime Exposure for the STS-ModAct trajectory -----#

        print('\nComputation of the Lifetime Exposure for the STS-ModAct trajectory')

        # if max threshold criteria met, run gmt mapping
        if d_isimip_meta[i]['GMT_STS_ModAct_valid']:

            print("ISMIP Simulation {} use for re-mapping GMT index for the STS-ModAct trajectory".format(i))

            # reindexing original exposure array based on GMT-mapping indices
            d_le_percountry_perrun_STS_ModAct = df_exposure_percountry.apply(
                lambda col: calc_life_exposure(
                    df_exposure_percountry.reindex(df_exposure_percountry.index[d_isimip_meta[i]['ind_RCP2GMT_STS_ModAct'][:]]).set_index(df_exposure_percountry.index),
                    df_life_expectancy_5,
                    col.name,
                ),
                axis=0,
            )

            # convert dataframe to data array of lifetime exposure (le) per country and birth year
            ds_le_percountry_perrun['le_percountry_perrun_STS_ModAct'].loc[{
                'run':i,
            }] = d_le_percountry_perrun_STS_ModAct.values.transpose() 

            #---------------------------------------------------------------------#
            # Per region                                                          #
            #---------------------------------------------------------------------#

            for region_ind, region in enumerate(ds_regions.region.values):
            
                region_name = ds_regions['name'].sel(region=region_ind).item()
                member_countries = ds_regions['member_countries'].sel(region=region_ind).values.tolist()

                weights = d_cohort_weights_regions[region_name]

                le_percountry_perrun = ds_le_percountry_perrun['le_percountry_perrun_STS_ModAct'].loc[{
                    'run': i,
                    'country': member_countries,
                    'birth_year': slice(None)
                }]

                weighted_avg = np.nansum(le_percountry_perrun.values * weights.values.T, axis=0) / np.nansum(weights.values.T, axis=0)

                ds_le_perregion_perrun['le_perregion_perrun_STS_ModAct'].loc[{
                    'run': i,
                    'region': region_ind,
                    'birth_year': slice(None)
                }] = weighted_avg
        
        else:

            print("ISMIP Simulation {} can not be used for re-mapping GMT for the STS-ModAct trajectory".format(i))

        #----- Computation of the Lifetime Exposure for the STS-Ren trajectory -----#

        print('\nComputation of the Lifetime Exposure for the STS-Ren trajectory')

        # if max threshold criteria met, run gmt mapping
        if d_isimip_meta[i]['GMT_STS_Ren_valid']:

            print("ISMIP Simulation {} use for re-mapping GMT index for the STS-Ren trajectory".format(i))

            # reindexing original exposure array based on GMT-mapping indices
            d_le_percountry_perrun_STS_Ren = df_exposure_percountry.apply(
                lambda col: calc_life_exposure(
                    df_exposure_percountry.reindex(df_exposure_percountry.index[d_isimip_meta[i]['ind_RCP2GMT_STS_Ren'][:]]).set_index(df_exposure_percountry.index),
                    df_life_expectancy_5,
                    col.name,
                ),
                axis=0,
            )

            # convert dataframe to data array of lifetime exposure (le) per country and birth year
            ds_le_percountry_perrun['le_percountry_perrun_STS_Ren'].loc[{
                'run':i,
            }] = d_le_percountry_perrun_STS_Ren.values.transpose() 

            #---------------------------------------------------------------------#
            # Per region                                                          #
            #---------------------------------------------------------------------#

            for region_ind, region in enumerate(ds_regions.region.values):
            
                region_name = ds_regions['name'].sel(region=region_ind).item()
                member_countries = ds_regions['member_countries'].sel(region=region_ind).values.tolist()

                weights = d_cohort_weights_regions[region_name]

                le_percountry_perrun = ds_le_percountry_perrun['le_percountry_perrun_STS_Ren'].loc[{
                    'run': i,
                    'country': member_countries,
                    'birth_year': slice(None)
                }]

                weighted_avg = np.nansum(le_percountry_perrun.values * weights.values.T, axis=0) / np.nansum(weights.values.T, axis=0)

                ds_le_perregion_perrun['le_perregion_perrun_STS_Ren'].loc[{
                    'run': i,
                    'region': region_ind,
                    'birth_year': slice(None)
                }] = weighted_avg
                    
        else:

            print("ISMIP Simulation {} can not be used for re-mapping GMT for the STS-Ren trajectory".format(i))

        #---------------------------------------------------------------------#
        # Computations of the Lifetime Exposure per country for the           #
        # Burning Embers diagram by mapping the GMTs of the ISIMIP            #
        # simulations to the stylized trajectories if the ISIMIP              #
        # simulations use is valid for remapping                              #
        #---------------------------------------------------------------------#   
        
        print('\nComputation of the Lifetime Exposure for the stylized trajectories of the BE')

        for step in GMT_labels:

            # if max threshold criteria met, run gmt mapping
            if d_isimip_meta[i]['GMT_strj_valid'][step]:

                print("ISMIP Simulation {} use for re-mapping GMT index = {} of the BE".format(i, step))
            
                # reindexing original exposure array based on GMT-mapping indices
                d_le_percountry_perrun_BE = df_exposure_percountry.apply(
                    lambda col: calc_life_exposure(
                        df_exposure_percountry.reindex(df_exposure_percountry.index[d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]).set_index(df_exposure_percountry.index),
                        df_life_expectancy_5,
                        col.name,
                    ),
                    axis=0,
                )
        
                # convert dataframe to data array of lifetime exposure (le) per country and birth year
                ds_le_percountry_perrun['le_percountry_perrun_BE'].loc[{
                    'run':i,
                    'GMT':step,
                }] = d_le_percountry_perrun_BE.values.transpose() 
    
                #---------------------------------------------------------------------#
                # Per region                                                          #
                #---------------------------------------------------------------------#

                for region_ind, region in enumerate(ds_regions.region.values):
                
                    region_name = ds_regions['name'].sel(region=region_ind).item()
                    member_countries = ds_regions['member_countries'].sel(region=region_ind).values.tolist()

                    weights = d_cohort_weights_regions[region_name]

                    le_percountry_perrun_BE = ds_le_percountry_perrun['le_percountry_perrun_BE'].loc[{
                        'run': i,
                        'GMT': step,
                        'country': member_countries,
                        'birth_year': slice(None)
                    }]

                    weighted_avg = np.nansum(le_percountry_perrun_BE.values * weights.values.T, axis=0) / np.nansum(weights.values.T, axis=0)

                    ds_le_perregion_perrun['le_perregion_perrun_BE'].loc[{
                        'run': i,
                        'GMT': step,
                        'region': region_ind,
                        'birth_year': slice(None)
                    }] = weighted_avg
            
            else:

                print("ISMIP Simulation {} can not be used for re-mapping GMT index = {} of the BE".format(i, step))

    # dump pickle of lifetime exposure per country
    with open(data_dir+'{}/{}/ds_le_percountry_perrun_gmt_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
        pk.dump(ds_le_percountry_perrun,f)
    
    # dump pickle of lifetime exposure per region
    with open(data_dir+'{}/{}/ds_le_perregion_perrun_gmt_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
        pk.dump(ds_le_perregion_perrun,f)

    return ds_le_percountry_perrun, ds_le_perregion_perrun
    
        
#%%---------------------------------------------------------------#
# Convert Area Fraction Affected (AFA) to                         #
# per-cohort number of extremes affecting one individual across   #
# life span to try to analyse the "age of emergence"              #
# by being time/age explicit to assess this. Did not come to      #
# fruition and not retain for further usage in Grant et al.(2025) #
#-----------------------------------------------------------------#

def calc_cohort_lifetime_exposure(
    d_isimip_meta,
    df_countries,
    countries_regions,
    countries_mask,
    da_population,
    da_cohort_size,
    flags,
):

    # loop over simulations
    for i in list(d_isimip_meta.keys()): 

        print('ISIMIP Simulation {} of {}'.format(i,len(d_isimip_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA = pk.load(f)

        # --------------------------------------------------------------------
        # per country 
        # --------------------------------------------------------------------

        # initialise dicts
        d_exposure_peryear_percountry = {}

        # get spatial average
        for j, country in enumerate(df_countries['name']):

            print('Processing country '+str(j+1)+' of '+str(len(df_countries))+' for Lifetime Exposure across cohorts', end='\r')
            
            #---------------------------------------------------#
            # calculate mean per country weighted by population #
            #---------------------------------------------------#

            # retrieve the indices of the countries following the list of df_countries in the countries_region regionmask.Regions
            # this indice is used after as an argument for the calc_weighted_fldmean() function. 
            ind_country = countries_regions.map_keys(country)

            # d_exposure_peryear_percountry is the dictionnary that for each country contains the data
            # array that contains for each time step (year) the exposure weighted by the population
            d_exposure_peryear_percountry[country] = calc_weighted_fldmean( 
                da_AFA,
                da_population, 
                countries_mask, 
                ind_country, 
                flag_region= False,
            )
            
            #---------------------------------------------------#

        # convert dictionary to data array
        da_exposure_peryear_percountry = xr.DataArray(
            list(d_exposure_peryear_percountry.values()),
            coords={
                'country': ('country', list(d_exposure_peryear_percountry.keys())),
                'time': ('time', da_AFA.time.values),
            },
            dims=[
                'country',
                'time',
            ],
        )

        # ------------------------------------------------------------------------------------#
        # v1 Not Used in Luke's framework : Compute da_exposure_cohort_strj by multiplying the 
        # da_exposure_peryear_percountry by da_cohort_size
        # -------------------------------------------------------------------------------------#

        # init of the DataArray that will contains the cohort exposure
        # for stylized trajectories

        # da_exposure_cohort_strj = xr.DataArray(
        #     coords={
        #         'country': ('country', list(d_exposure_peryear_percountry.keys())),
        #         'time': ('time', da_AFA.time.values),
        #         'ages': ('ages', da_cohort_size.ages.values),
        #         'GMT': ('GMT', GMT_labels),
        #     },
        #     dims=[
        #         'country',
        #         'time',
        #         'ages',
        #         'GMT',
        #     ],
        # ) 
        
        # # GMT mapping for cohort exposure for stylized trajectories
        # # with corresponding historical + RCP simulations
        # for step in GMT_labels:

        #     # if max threshold criteria met, run gmt mapping              
        #     if d_isimip_meta[i]['GMT_strj_valid'][step]:
                
        #         # assign data to data array based on step in stylized trajectories
        #         da_exposure_cohort_strj.loc[
        #             {'country':da_exposure_cohort_strj.country,
        #                 'time':da_exposure_cohort_strj.time,
        #                 'GMT':step}
        #             ] = da_exposure_peryear_percountry.reindex(
        #                 {'time':da_exposure_peryear_percountry['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
        #             ).assign_coords({'time':year_range}) * da_cohort_size

        # with open(data_dir+'{}/{}/exposure_cohort_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],i), 'wb') as f:
        #     pk.dump(da_exposure_cohort_strj,f)
            
        #-------------------------------------------------------------------------------#
        # v2 Used in Luke's framework : Compute da_exposure_peryear_perage_percountry_strj 
        # by multiplying the da_exposure_peryear_percountry 
        # by xr.full_like(da_cohort_size,1). Compared to v1 it includes an 'ages' 
        # dimension which makes the v2 more detailed. 
        #-------------------------------------------------------------------------------#

        # init of the DataArray that will contains the cohort exposure
        # for stylized trajectories
        da_exposure_peryear_perage_percountry_strj = xr.DataArray(
            coords={
                'country': ('country', list(d_exposure_peryear_percountry.keys())),
                'time': ('time', da_AFA.time.values),
                'ages': ('ages', da_cohort_size.ages.values),
                'GMT': ('GMT', GMT_labels),
            },
            dims=[
                'country',
                'time',
                'ages',
                'GMT',
            ],
        )
        
        #------------------------------------------------------------------------------------#
        # GMT mapping with corresponding historical + RCP simulations 
        # for stylized trajectories in dimension expansion of da_exposure_peryear_percountry
        
        for step in GMT_labels:
            
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                
                # assign data to data array based on step in stylized trajectories
                da_exposure_peryear_perage_percountry_strj.loc[
                    {'country':da_exposure_peryear_perage_percountry_strj.country,
                        'time':da_exposure_peryear_perage_percountry_strj.time,
                        'GMT':step}
                    ] = da_exposure_peryear_percountry.reindex(
                        {'time':da_exposure_peryear_percountry['time'][d_isimip_meta[i]['ind_RCP2GMT_strj'][:,step]]}
                    ).assign_coords({'time':year_range}) * xr.full_like(da_cohort_size,1)
        
        #------------------------------------------------------------------------------------#

        with open(data_dir+'{}/{}/exposure_peryear_perage_percountry_{}_{}.pkl'.format(flags['version'],flags['extr'],flags['extr'],i), 'wb') as f:
            pk.dump(da_exposure_peryear_perage_percountry_strj,f)                     
        
#%%---------------------------------------------------------------#
# Convert PIC Area Fraction Affected (AFA) to                     #
# per-country number of extremes affecting one individual         #
# across life span under PIC conditions for climate               #
# and 1960 demography                                             #
#-----------------------------------------------------------------#

def calc_lifetime_exposure_pic(
    d_pic_meta, 
    df_countries, 
    countries_regions, 
    countries_mask, 
    da_population, 
    df_life_expectancy_5, 
    flags,
):

    d_pic_le_percountry_perrun = {}                 
    
    # loop over simulations
    for n,i in enumerate(list(d_pic_meta.keys())):

        print('PIC Simulation '+str(n+1)+ ' of '+str(len(d_pic_meta)))

        # load AFA data of that run
        with open(data_dir+'{}/{}/isimip_AFA_pic_{}_{}.pkl'
        .format(flags['version'],flags['extr'],flags['extr'],str(i)), 'rb') as f:
            da_AFA_pic = pk.load(f)
        
        # get 1960 life expectancy in each country
        life_expectancy_1960 = xr.DataArray(
            df_life_expectancy_5.loc[1960].values,
            coords={
                'country': ('country', df_life_expectancy_5.columns)
            }
        )       
        
        # --------------------------------------------------------------------
        # per country 
        # --------------------------------------------------------------------
        
        # initialise dicts
        d_exposure_peryear_percountry_pic = {}
        
        # get spatial average
        for j, country in enumerate(df_countries['name']): # Luke's comment : with other stuff running, this loop took 91 minutes therefore consider first doing the weighted mean and then boot strapping? does that make sense?

            print('Processing country '+str(j+1)+' of '+str(len(df_countries))+' for Lifetime Exposure across cohorts', end='\r')

            #---------------------------------------------------#
            # calculate mean per country weighted by population #
            #---------------------------------------------------#

            # retrieve the indices of the countries following the list of df_countries in the countries_region regionmask.Regions
            # this indice is used after as an argument for the calc_weighted_fldmean() function.
            ind_country = countries_regions.map_keys(country)

            
            # d_exposure_peryear_percountry_pic is the dictionnary that for each country contains the data
            # array that contains for each time step (year) the exposure weighted by the population
            # perform for corresponding picontrol climate conditions 
            # and assuming constant 1960 population density (this line takes about 16h by itself)
            d_exposure_peryear_percountry_pic[country] = calc_weighted_fldmean(
                da_AFA_pic, 
                da_population[0,:,:], # use of the earliest year used for population weights
                countries_mask, 
                ind_country, 
                flag_region= False,
            )
            
            #---------------------------------------------------#

        # convert dictionary to data array
        da_exposure_pic = xr.DataArray(
            list(d_exposure_peryear_percountry_pic.values()),
            coords={
                'country': ('country', list(d_exposure_peryear_percountry_pic.keys())),
                'time': ('time', da_AFA_pic.time.values),
            },
            dims=[
                'country',
                'time',
            ],
        )

        # bootstrap native pic exposed area data ;pic_life_extent, nboots, resample_dim
        da_exposure_pic = xr.concat([resample(da_exposure_pic,resample_dim,pic_life_extent) for i in range(nboots)],dim='lifetimes')

        # ------------------------------------------------------------------------- #
        # calculate the cumulative exposure for lifetime under pic climate conditions
        # for the 1960 demographics in all countries
        # ------------------------------------------------------------------------- #
        # add the exposure for all years that will be lived by the individuals of a 
        # particular cohort with regards to their birth year and their associated 
        # life expectancy. Add the fraction of the exposure of their last year of life
        d_pic_le_percountry_perrun[i] = da_exposure_pic.where(da_exposure_pic.time < 1960 + np.floor(life_expectancy_1960)).sum(dim='time') + \
            da_exposure_pic.where(da_exposure_pic.time == 1960 + np.floor(life_expectancy_1960)).sum(dim='time') * \
                (life_expectancy_1960 - np.floor(life_expectancy_1960))
        # ------------------------------------------------------------------------- #

    # save pickles
    with open(data_dir+'{}/{}/d_pic_le_percountry_perrun.pkl'.format(flags['version'],flags['extr']), 'wb') as f:
        pk.dump(d_pic_le_percountry_perrun,f)

    return d_pic_le_percountry_perrun

#%%---------------------------------------------------------------#
# Computes the Exposure Multiplication Factor (EMF)               #
# based on ms_exposure.m from Thiery et al.(2021)                 #
# Take as input the dataset containing the MMM lifetime exposure  #
# Option to use the MMM on the historical period for the          #
# reference level of exposure or the MMM of the picontrol         #
# simulations                                                     #
#-----------------------------------------------------------------#

def calc_EMF(
    flags,
    ds_le_exposure, 
    ref_pic,
):

    #--------------------- Definition of the reference exposure -----------------#

    # Decide which reference to use: multi-model mean of the picontrol simulations (computed in ???) 
    # or of the historical simulations (computed in this function)
    if ref_pic:
        exposure_ref = exposure_pic_mean
        
    else:
        exposure_ref = ds_le_exposure

    # Step 1: extract the 'mmm' variable
    da_ref = exposure_ref['mmm_BE']

    # Step 2: select only the reference year along the 'birth_year' dimension
    da_ref_slice = da_ref.sel(birth_year=year_ref)

    # Step 3: broadcast the reference slice across all birth years
    birth_years = ds_le_exposure.coords['birth_year']
    da_ref_broadcasted = da_ref_slice.expand_dims(birth_year=birth_years)

    #----------------- Retrieval of the lifetime exposure values -----------------#

    if 'region' in ds_le_exposure.dims:

        da_mmm = ds_le_exposure['mmm_BE'].sel(
        birth_year=ds_le_exposure.coords['birth_year'],
        region=ds_le_exposure.coords['region'],
        GMT=ds_le_exposure.coords['GMT']
        )

    if 'country' in ds_le_exposure.dims:

        da_mmm = ds_le_exposure['mmm_BE'].sel(
        birth_year=ds_le_exposure.coords['birth_year'],
        country=ds_le_exposure.coords['country'],
        GMT=ds_le_exposure.coords['GMT']
        )

    # Option : drop any remaining dimensions (like 'run') if needed
    da_mmm = da_mmm.squeeze(drop=True)    

    #----------------------------- EMF Computation --------------------------------#

    # At this point, da_ref_broadcasted has the same shape as da_mmm

    ds_EMF_mmm = da_mmm / da_ref_broadcasted

    # Set EMF of infinitue value to 100 #

    # Calculate the percentage of infinite values
    ninf = np.isinf(ds_EMF_mmm).sum().item() / ds_EMF_mmm.size * 100  # e.g., 2.8%

    # Replace infinite values by 100
    ds_EMF_mmm = ds_EMF_mmm.where(~np.isinf(ds_EMF_mmm), 100)

    # Save pickles
    if 'region' in ds_le_exposure.dims:
        with open(data_dir+'{}/{}/ds_EMF_perregion_gmt_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
            pk.dump(ds_EMF_mmm,f)
    
    if 'country' in ds_le_exposure.dims:
        with open(data_dir+'{}/{}/ds_EMF_percountry_gmt_{}.pkl'.format(flags['version'],flags['extr'],flags['gmt']), 'wb') as f:
            pk.dump(ds_EMF_mmm,f)

    print("\nEnd of the computation\n")

    return ds_EMF_mmm