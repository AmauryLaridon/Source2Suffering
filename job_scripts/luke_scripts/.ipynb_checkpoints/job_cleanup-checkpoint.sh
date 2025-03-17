#!/bin/bash

#SBATCH --job-name=cleanup
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.grant@vub.be
#SBATCH --output out_cleanup
#SBATCH --error err_cleaup

echo Start Job
date

cd /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00012/vsc10116/lifetime_exposure_isimip
./cleanup.sh

echo End Job
