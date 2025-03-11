#!/bin/bash

#SBATCH --job-name=cropfailedarea
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.grant@vub.be
#SBATCH --output out_cropfailedarea
#SBATCH --error err_cropfailedarea

echo Start Job
date

cd /data/brussel/vo/000/bvo00012/vsc10116/lifetime_exposure_isimip

# load modules
# ml Python/3.9.6-GCCcore-11.2.0
module load Python/3.10.4-GCCcore-11.3.0 # new

# ml geopandas/0.11.0-foss-2021b
# module load geopandas/0.14.2-foss-2023a
module load geopandas/0.12.2-foss-2022a

# ml openpyxl/3.0.9-GCCcore-11.2.0
module load openpyxl/3.0.10-GCCcore-11.3.0

# ml regionmask/0.9.0-foss-2021b
module load regionmask/0.9.0-foss-2022a

# ml xarray/0.20.1-foss-2021b

module load xarray/2022.6.0-foss-2022a # can keep as is for now

# ml netcdf4-python/1.5.8-foss-2021b
# module load netcdf4-python/1.6.4-foss-2023a
module load netcdf4-python/1.6.1-foss-2022a

# ml SciPy-bundle/2021.10-foss-2021b # keep as is for now
module load SciPy-bundle/2022.05-foss-2022a

# ml matplotlib/3.4.3-foss-2021b
# module load matplotlib/3.7.2-gfbf-2023a
module load matplotlib/3.5.2-foss-2022a

# ml Cartopy/0.20.3-foss-2021b
module load Cartopy/0.20.3-foss-2022a

# ml Shapely/1.8.2-foss-2021b
module load Shapely/1.8.2-foss-2022a
# module swap xarray/0.20.1-foss-2021b xarray/2022.6.0-foss-2021b

python3 main.py

echo End Job
