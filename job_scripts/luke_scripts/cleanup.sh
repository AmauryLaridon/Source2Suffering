#!/bin/bash

pickle_dir=/vscmnt/brussel_pixiu_data/_data_brussel/vo/000/bvo00012/vsc10116/lifetime_exposure_isimip/data/pickles
cd $pickle_dir

extrs_list=("heatwavedarea" "cropfailedarea" "driedarea" "floodedarea" "tropicalcyclonedarea" "burntarea")

# all pi files for this mod
c_files=($(find ./ -maxdepth 1 -mindepth 1 -type f -name "gridscale_dmg_*" -printf '%P\n'))
countries=()

# retrieve realisation (member) labels
for cf in "${c_files[@]}"; do

    c_pkl="$(cut -d '_' -f 3 <<<"$cf")"
    c="$(cut -d'.' -f1 <<<"$c_pkl")"
    countries[${#countries[@]}]=$c
    
done

for extr in "${extrs_list[@]}"; do

    cd $extr
    echo "removing intermediary files for ${extr}"

    rm exposure_cohort*
    rm exposure_peryear_perage_percountry*
    rm ds_exposure_aligned_"${extr}"*
    rm da_emergence_mask_"${extr}"*
    rm ds_age_emergence_"${extr}"*
    rm da_age_emergence_"${extr}"*

    for c in "${countries[@]}"; do

        rm gridscale_le_pic_"${extr}"_"${c}"*
        rm gridscale_le_"${extr}"_"${c}"*
        rm gridscale_emergence_mask_"${extr}"_"${c}"*
        rm gridscale_spatially_explicit_"${extr}"_"${c}"*
        rm gridscale_age_emergence_"${extr}"_"${c}"*

    done

    cd $pickle_dir

done

