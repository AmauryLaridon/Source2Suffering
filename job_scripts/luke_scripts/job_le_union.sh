#!/bin/bash

#SBATCH --job-name=union
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.grant@vub.be
#SBATCH --output out_union
#SBATCH --error err_union

echo Start Job
date

cd /data/brussel/vo/000/bvo00012/vsc10116/lifetime_exposure_isimip

# load modules
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

python3 main.py

echo End Job
