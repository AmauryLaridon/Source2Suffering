#!/bin/bash

#SBATCH --job-name=heatwavedarea
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.grant@vub.be
#SBATCH --output out_heatwavedarea
#SBATCH --error err_heatwavedarea

echo Start Job
date

cd /data/brussel/vo/000/bvo00012/vsc10116/lifetime_exposure_isimip

# load modules
module load Python/3.9.6-GCCcore-11.2.0
module load geopandas/0.11.0-foss-2021b
module load openpyxl/3.0.9-GCCcore-11.2.0
module load regionmask/0.9.0-foss-2021b
#ml xarray/0.20.1-foss-2021b
module load xarray/2022.6.0-foss-2022a
module load netcdf4-python/1.5.8-foss-2021b
module load SciPy-bundle/2021.10-foss-2021b
module load matplotlib/3.4.3-foss-2021b
module load Cartopy/0.20.3-foss-2021b
module load Shapely/1.8.2-foss-2021b
#module swap xarray/0.20.1-foss-2021b xarray/2022.6.0-foss-2021b

python3 main.py

echo End Job
