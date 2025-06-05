#!/bin/bash

#SBATCH --job-name=s2s_tropicalcyclonedarea
#SBATCH --ntasks=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amaury.laridon@vub.be
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

timestamp=$(date +"%Y%m%d_%H%M")

output_file="/data/brussel/vo/000/bvo00012/vsc11137/source2suffering/output/job_output/out_s2s_tropicalcyclonedarea_${timestamp}"
error_file="/data/brussel/vo/000/bvo00012/vsc11137/source2suffering/output/job_error/err_s2s_tropicalcyclonedarea_${timestamp}"

exec > >(tee -a "$output_file") 2> >(tee -a "$error_file" >&2)

echo -------------------------------------------------------------------------------------------------

echo "Start Job - Job Name: $SLURM_JOB_NAME - Job ID: $SLURM_JOB_ID"

date

echo "--- System Hardware Info ---"

echo "Running on partition: $SLURM_JOB_PARTITION"

echo -n "CPU Model: "; lscpu | grep "Model name" | awk -F': ' '{print $2}'

echo -n "Number of CPU Cores: "; nproc

echo -n "Total RAM: "; free -h | awk '/^Mem:/ {print $2}'

echo -n "Disk Space: "; df -h --total | grep "total" | awk '{print $2}'

echo -n "GPU Info: "; lspci | grep -E "VGA|3D" | awk -F': ' '{print $2}'

echo -n "System Architecture: "; uname -m


echo -------------------------------------------------------------------------------------------------


cd /data/brussel/vo/000/bvo00012/vsc11137/source2suffering

echo Start loading the modules for the Source2Suffering Project

#### load modules ###

start_time_mod=$(date +%s.%N)  # capture the start time

module load Python/3.10.4-GCCcore-11.3.0 
module load geopandas/0.12.2-foss-2022a
module load openpyxl/3.0.10-GCCcore-11.3.0
module load regionmask/0.9.0-foss-2022a
module load xarray/2022.6.0-foss-2022a 
module load netcdf4-python/1.6.1-foss-2022a
module load SciPy-bundle/2022.05-foss-2022a
module load matplotlib/3.5.2-foss-2022a
module load Cartopy/0.20.3-foss-2022a
module load Shapely/1.8.2-foss-2022a
module load Seaborn/0.12.1-foss-2022a
module load imageio/2.22.2-foss-2022a 
module load ncview/2.1.8-gompi-2022a

end_time_mod=$(date +%s.%N)  # Capture the end time
execution_time_mod=$(echo "$end_time_mod - $start_time_mod" | bc)
execution_time_mod_min=$(echo "scale=2; $execution_time_mod / 60" | bc)

echo "Execution time to load the modules: $execution_time_mod_min minutes"

echo -------------------------------------------------------------------------------------------------

#### Execute Scripts ###

echo Start executing the scripts

start_time_main=$(date +%s.%N)

export EXTR_VALUE="tropicalcyclonedarea"

echo "Running main.py with extr=$EXTR_VALUE"

# Run the main script without redirection (it's already done via exec)
python3 main.py

end_time_main=$(date +%s.%N)
execution_time_main=$(echo "$end_time_main - $start_time_main" | bc)
execution_time_main_min=$(echo "scale=2; $execution_time_main / 60" | bc)

echo End Job 
date
echo -------------------------------------------------------------------------------------------------
