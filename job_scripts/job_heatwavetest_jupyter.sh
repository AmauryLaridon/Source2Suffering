#!/bin/bash

#SBATCH --job-name jupyter 
#SBATCH --output jupyter.log 
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.grant@vub.be

echo $HOSTNAME
export OMP_NUM_THREADS=1

# load modules
ml JupyterLab/3.1.6-GCCcore-11.2.0
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
ml dask/2022.1.0-foss-2021b
module swap xarray/0.20.1-foss-2021b xarray/2022.6.0-foss-2021b
jupyter-lab --no-browser 

