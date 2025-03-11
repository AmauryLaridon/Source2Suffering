#!/bin/bash

#SBATCH --job-name=copydata
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amaury.laridon@vub.be
#SBATCH --output out_copydata
#SBATCH --error err_copydata

echo Start Job

date

rsync -avh --progress --partial --append /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00012/vsc10116/lifetime_exposure_isimip/data/ /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00012/vsc11137/source2suffering/data

echo End Job
