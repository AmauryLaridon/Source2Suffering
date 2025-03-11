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
# based on Luke job script job_le.sh

ml Python/3.9.6-GCCcore-11.2.0
ml geopandas/0.11.0-foss-2021b
ml openpyxl/3.0.9-GCCcore-11.2.0
ml regionmask/0.9.0-foss-2021b
ml xarray/0.20.1-foss-2021b
ml netcdf4-python/1.5.8-foss-2021b
ml SciPy-bundle/2021.10-foss-2021b
ml matplotlib/3.4.3-foss-2021b
ml Cartopy/0.20.3-foss-2021b
ml Shapely/1.8.2-foss-2021b
module swap xarray/0.20.1-foss-2021b xarray/2022.6.0-foss-2021b

date
echo End Job 
