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
# based on Luke job script job_le_burntarea_newmodules.sh

ml Python/3.10.4-GCCcore-11.3.0 geopandas/0.12.2-foss-2022a openpyxl/3.0.10-GCCcore-11.3.0 regionmask/0.9.0-foss-2022a xarray/2022.6.0-foss-2022a netcdf4-python/1.6.1-foss-2022a SciPy-bundle/2022.05-foss-2022a matplotlib/3.5.2-foss-2022a Cartopy/0.20.3-foss-2022a Shapely/1.8.2-foss-2022a

date
echo End Job 
