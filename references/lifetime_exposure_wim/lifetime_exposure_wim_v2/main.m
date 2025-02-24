

% --------------------------------------------------------------------
% main script to postprocess and visualise ISIMIP2b output: 'Lange data'
% For 'burden of young people' study
% --------------------------------------------------------------------


% to do:
% - is averaging across models and extremes ok for EMF??? Improve uncertainty propagation?
% - use spatially-explicit life expectancy data for BE too? (now it's using the global LE)
% - avoid double counting of ISIMIP years ?
% - make maps of absolute changes
% - regional line plots: add pie charts with relative fractions for selected age groups (at birth? end of life is more difficult and SSP dependend) below the figure
% - ms_exposure: mf_invprctile and prctile not fully reversible !!!
% - BE scenarios branch off in 2010 already, is this ok? We've had a decade of warming since then already, so most conservative scenario is not possible anymore
% - apply running mean to GMT data? e.g. 5-yr running mean? to be able to reduce maxdiff from 0.5 to 0.1 again?


% - wildfires, global average is higher than that in any region!!!! can this be?
% - 405 vs 393 nruns, see input to mf_exposure_mmm function (TCs seem to be copied twice! but I don't think in this case the mmm is adversely affected!)


% notes:
% - ISIMIP GMT anomalies are computed using the 1850-1900 period (51-yr average) as reference (historical, not picontrol), consistent with SR15 Ch1 (sect. 1.2.1)
% - in several countries life expectancy in 1960 is <60, yet x-axis of contry plots suggests they are still alive in 2020!
% - previously, hazards were dropping for very young ages. This was due to makima giving weird extrapolations in some countries. Linear extrapolation solved the issue
% - hazards emerge in certain countries, e.g. TCs in Korea,
% - heatwaves suspicious in eastern-Europe+scandinavia+UK: seems to be related to fact that definition is for 'wet heatwaves' see e-mail stefan Lange 20/01/2020
% - population data assumes SSP3 for all RCP/GMT scenarios
% - spatial averaging is done using world life expectancy instead of country life expectancy
% - now brute GMT difference threshold (0.5K) to remove entire runs, this could be refined if needed
% - updated axdiff threshold from 0.1°C to 0.5°C due to presence of year-to-year jumps >0.1°C in annual GMT series (and sticking to annual
%   makes sense due to damage function approach + moving window has no data for last years)
% - for the global analysis we assume one uniform life expectancy, but instead we could account for the spatial life expectancy variability
% - BE: weird line structures are due to GCM sampling.
% - spatial pdfs: mismatch of ~500 million people when applying repelem at pixel scale using country-level data because they live along the coastlines
%     popmask = zeros(size(population(:,:,61)));
%     for i=1:ncountries; popmask(countries.mask{i})=1; end;
%     figure;imagesc(popmask);colorbar
%     AAE=population(:,:,61);
%     figure;imagesc(AAE);colorbar;caxis([0 1]);
% - stefan's new data set (20200417) differs from his earlier data set (20190601). We can use new data because:
%     He uses the same heatwave definition but processed it in cdo instead of R to be faster
%     AFA=NaN instead of AFA=zero where there is no data ==> Small islands like Bahama's, Comores, Cape Verde etc. now sometimes have NaN as e.g. CLM4.5, LPJmL,... do not simulate there. So caution needed when interpreting small island results as there is less data
%     the original files differ, e.g. burntarea orchidee-gfdl-rcp2.6. So it looks like the differences are owing to new input data, not flaws in my procesisng chain
% - regarding means:
%     geometric mean on exposure suggests no change in many countries: does not make any sense if you look at EMFs of individual hazards ==> don't use
%     geometric mean on EMFs : works, gives substantial diff between 1.5 and NDC, is theoretically defendable (mean of 'different things')
%     harmonic mean on EMFs : works and is theoretically defendable (mean of ratios), but is most conservative of all means and gives little diff between scenarios and weird results in many countries under 1.5



tic


% clean up
clc;
clear;
close all;


% flags
flags.extr  = 0;        % 0: all
                        % 1: burntarea
                        % 2: cropfailedarea
                        % 3: driedarea
                        % 4: floodedarea
                        % 5: heatwavedarea
                        % 6: tropicalcyclonedarea
flags.masks = 0;        % 0: do not process country data (i.e. load masks workspace)
                        % 1: process country data (i.e. produce and save masks as workspace)
flags.runs  = 1;        % 0: do not process ISIMIP runs (i.e. load runs workspace)
                        % 1: process ISIMIP runs (i.e. produce and save runs as workspace)
flags.exposure = 1;     % 0: do not process ISIMIP runs to compute exposure (i.e. load exposure workspace)
                        % 1: process ISIMIP runs to compute exposure (i.e. produce and save exposure as workspace)
flags.exposure_pic = 1; % 0: do not process ISIMIP runs to compute picontrol exposure (i.e. load exposure workspace)
                        % 1: process ISIMIP runs to compute picontrol exposure (i.e. produce and save exposure as workspace)
flags.coldwaves = 0;    % 0: compute heatwavedarea (default)
                        % 1: compute coldwavedarea instead of heatwavedarea
flags.plot  = 1;        % 0: do not plot
                        % 1: plot
flags.valp  = 1;        % 0: do not compute values used in the paper
                        % 1: compute values used in the paper
flags.valc  = 1;        % 0: do not compute values used for communication
                        % 1: compute values used for communication
flags.saveall = 1;      % 0: do not save all data in matlab workspace
                        % 1: save all data in matlab workspace



% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% declare globals
global island                                                              %#ok<NUSED>


% add matlab scripts directory to path
%addpath(genpath('C:\Users\Wim Thiery\Documents\research\matlab_scripts'));
addpath(genpath('C:\Users\AL000111\Documents\Source2Suffering\references\lifetime_exposure_wim\matlab_scripts'));


% add directory containing nc files to path
% indir = 'C:\Users\Wim Thiery\Documents\research\ISIMIP2b_exposure\ncfiles';
indir = 'C:\Users\AL000111\Documents\Source2Suffering\references\lifetime_exposure_wim\lifetime_exposure_wim_v2\ncfiles';
addpath(genpath(indir));  


% initialise age and associated time period of interest
ages        = (60:-1:0)';
age_young   = 0;
age_ref     = max(ages);
year_ref    = 2020;
year_start  = year_ref - age_ref;
birth_years = (year_start:year_ref)';       
year_end    = 2113;                   % based on maximum life expectancy reported in UN WPP


% initialise age groups
% (https://www.carbonbrief.org/analysis-why-children-must-emit-eight-times-less-co2-than-their-grandparents)
% (https://www.pewresearch.org/fact-tank/2019/01/17/where-millennials-end-and-generation-z-begins/)
agegroups = {'Boomers'    1950 1965;
             'Gen X'      1965 1981;
             'Millenials' 1981 1997;
             'Gen Z'      1997 2020;};

         
% initialise reference period for computing GMT anomalies
year_start_GMT_ref  = 1850;
year_end_GMT_ref    = 1900;


% initialise types of extremes
extremes        = {'burntarea', 'cropfailedarea', 'driedarea', 'floodedarea' , 'heatwavedarea', 'tropicalcyclonedarea'};
extremes_legend = {'Wildfires', 'Crop failures' , 'Droughts' , 'River floods', 'Heatwaves'    , 'Tropical cyclones'   , 'All'};
if flags.extr > 0
    extremes        = extremes(flags.extr);
    extremes_legend = extremes_legend(flags.extr);
end


% initialise model names
model_names.burntarea            = {'CARAIB', 'LPJ-GUESS', 'LPJmL', 'ORCHIDEE', 'VISIT'                                        };
model_names.cropfailedarea       = {'GEPIC' , 'LPJmL'    , 'PEPIC'                                                             };
model_names.driedarea            = {'CLM45' , 'H08'      , 'LPJmL', 'JULES-W1', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'WaterGAP2'};
model_names.floodedarea          = {'CLM45' , 'H08'      , 'LPJmL', 'JULES-W1', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'WaterGAP2'};
model_names.heatwavedarea        = {'HWMId-humidex', 'HWMId99-humidex40', 'HWMId97p5-humidex40', 'HWMId99-tasmax35', 'HWMId97p5-tasmax35', 'HWMId99', 'HWMId97p5', 'humidex40d3', 'humidex40d5', 'humidex45d3', 'humidex45d5', 'CWMId99'};
model_names.tropicalcyclonedarea = {'KE-TG-meanfield'};


% Set threshold maximum T difference between RCP and GMT trajectories
% i.e. any run with T difference exceeding this threshold is excluded
% year-to-year jumps in GMT larger than 0.1, so using a 0.1 maxdiff threshold erronously removes runs
% used to be 0.5, but then rcp2.6 is used for high-warming levels
% Anything between 0.1 and 0.2 removes RCP2.6 in NDC scenarios (see histograms of maxdiff_NDC)
% take 0.2 to have more data in BE scenarios and hence smooth EMF curves in BE plot
RCP2GMT_maxdiff_threshold = 0.2; % [K]


% set kernel x-values
kernel_x = 1:0.5:50;


% Define Transient Climate response to cumulative emissions.
% Expert advice (27/03/2024): Je kan 1.65°C per 1000 PgC gebruiken, maar (mocht dat nuttig zijn) kan je ook transparant gecommuniceerde hogere percentielen gebruiken. 
% Die geven namelijk een risicoperspectief. Bijvoorbeeld, 2.2°C per 1000 PgC is nog steeds enkel het 66ste percentiel. Dat is niet per se extreem.
% Deze methode is wel vooral toepasbaar op CO2 and niet CO2-eq. Dus als er veel methaan in die emissies zit dan zou je dat een beetje moeten aanpassen.
TCRE = 0.45 ./ 1000E9; % 1.65°C / 3.7 = 0,45°C per 1000 Gt CO2eq


% define the mortality cost of Carbon
mortality_cost_carbon = 4434; % “4,434 metric tons of carbon dioxide in 2020 […] causes one excess death globally in expectation between 2020-2100” (Bressler, 2021)



% --------------------------------------------------------------------
% load data
% --------------------------------------------------------------------


ms_load



% --------------------------------------------------------------------
% manipulations: general
% --------------------------------------------------------------------


ms_manip



% --------------------------------------------------------------------
% Prepare for burning embers figure
% --------------------------------------------------------------------


ms_burning_embers;



% --------------------------------------------------------------------
% Compute exposure per lifetime
% --------------------------------------------------------------------


ms_exposure;



% --------------------------------------------------------------------
% visualise output
% --------------------------------------------------------------------


if flags.plot == 1
   ms_plotscript
end



% --------------------------------------------------------------------
% get values used in the paper
% --------------------------------------------------------------------


if flags.valp == 1
   ms_valp
end



% --------------------------------------------------------------------
% get values used for communication
% --------------------------------------------------------------------


if flags.valc == 1
   ms_valc
end



% --------------------------------------------------------------------
% get values used in STC report
% --------------------------------------------------------------------


if flags.saveall == 1
    
    % save arrays as matlab workspace
    disp('saving mw_output')
    save('mw_output');

end


toc
