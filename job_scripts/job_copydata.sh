#!/bin/bash

$SBATCH --job-name=copydata
$SBATCH --time=24:00:00
$SBATCH --ntasks=1
$SBATCH --mem=32G
$SBATCH --mail-type=ALL
$SBATCH --mail-user=amaury.laridon@vub.be
$SBATCH --output out_copydata
$SBATCH --error err_copydata

echo Start Job

date

rsync -avh --progress --partial --append /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00012/vsc11137/source2suffering/references/lifetime_exposure_luke/lifetime_exposure_luke_original/data/ /vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00012/vsc11137/source2suffering/data

echo End Job
