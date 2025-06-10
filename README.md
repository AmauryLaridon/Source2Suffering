# Source2Suffering

Repository containing the scripts and data needed for the Source2Suffering project. The scripts and data are based on L.Grant et al.(2025) and Thiery et al.(2021) 

## How To Run ?

### Running from VUB On Demand for HYDRA HPC

1. Before launching the VS Code server you have to paste the following command on the "Pre-run sciptlet" :

ml Python/3.10.4-GCCcore-11.3.0 geopandas/0.12.2-foss-2022a openpyxl/3.0.10-GCCcore-11.3.0 regionmask/0.9.0-foss-2022a xarray/2022.6.0-foss-2022a netcdf4-python/1.6.1-foss-2022a SciPy-bundle/2022.05-foss-2022a matplotlib/3.5.2-foss-2022a Cartopy/0.20.3-foss-2022a Shapely/1.8.2-foss-2022a Seaborn/0.12.1-foss-2022a

2. Once the VS Code server has been launched, select the corresponding Python 3.10.4 64-bit interpreter
3. Open the _settings.py_ script and change the _scripts_dir_ variable with the exact location of your scripts on Hydra.
4. To avoid a _CryptographyDeprecationWarning: Blowfish has been deprecated_ warning in the output you can execute the following command : _pip install --upgrade paramiko cryptography_
5. Configure all flags based on the analysis you would like to perform. The _settings.py_ scripts offer other configuration details that can be change.
6. Execute the main.py script.

## Version used for the Assessments 

 **Born into the Climate Crisis 2 (BiCC2) - Norway** : Commit full SHA : 407ee6bb38a0d1bbcf4918837fa8d675929431fd

 **Born into the Climate Crisis 2 (BiCC2) - Spain** : Commit full SHA : 37864688cbe3b507bbcdd258f9cd2c9e16dbd8b3

 **Nikkei Journal - Japan** : Commit Full SHA : 
 7b4f99ce4c5b150d82d3f5a773c83d588ff21f75
