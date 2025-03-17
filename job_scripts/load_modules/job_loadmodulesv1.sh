#!/bin/bash

#SBATCH --job-name=load_modules 
#SBATCH --output /data/brussel/vo/000/bvo00012/vsc11137/source2suffering/output/out_loadmodules
#SBATCH --error /data/brussel/vo/000/bvo00012/vsc11137/source2suffering/output/error/err_load_modules
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amaury.laridon@vub.be

echo Start loading the modules for the Source2Suffering Project
date

# load modules
# based on the latest versions of the modules available on Hydra

ml xarray/2023.9.0-gfbf-2023a
ml Cartopy/0.22.0-foss-2023a
ml openpyxl/3.1.2-GCCcore-12.3.0
ml geopandas/0.14.2-foss-2023a
ml netcdf4-python/1.6.4-foss-2023a
ml regionmask/0.12.1-foss-2023a
ml Dask-ML/2024.4.4-foss-2023a
ml Shapely/2.0.1-gfbf-2023a

date
echo End Job 
