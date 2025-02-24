

% --------------------------------------------------------------------
% DEM of selected domain
% --------------------------------------------------------------------


tic


% clean up
clc;
clear;
close all;
          
               
% flags
flag_EAfr = 1; % 0: do not plot SRTM DEM
               % 1: plot SRTM-DEM
flag_LVB  = 0; % 0: do not plot big domain
               % 1: plot big domain

               
               
% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% add matlab scripts directory to path
addpath(genpath('C:\Users\u0079068\Documents\Research\matlab_scripts'));



% --------------------------------------------------------------------
% DEM of selected domain
% --------------------------------------------------------------------


if flag_EAfr == 1

    
% definegrid boundaries
lat_min_EAfr = -17; % = -3.5 - 6.875 - 5;
lat_max_EAfr = 11;  % = -3.6 + 6.875 + 5;
lon_min_EAfr = 15;  % = 31 - 5.625 - 5;
lon_max_EAfr = 51;  % = 31 + 5.625 + 5;


% call import and plot function
X = mf_readhgt(lat_min_EAfr:lat_max_EAfr, lon_min_EAfr:lon_max_EAfr, 'merge', 'plot');


% draw 12 km domain borders
lat_min    = -12; % 163 * 0.11°
lat_max    = 6;
lon_min    = 22;  % 200 * 0.11°
lon_max    = 44;
rectangle('Position',[lon_min   , lat_min   , lon_max    - lon_min   , lat_max    - lat_min   ],'LineWidth',2.5)


% draw 2.8 km domain borders
lat_min_in = -5.0; % 200 * 0.025°
lat_max_in =  2.5;
lon_min_in = 27.5; % 440 * 0.025°
lon_max_in = 38.5;
rectangle('Position',[lon_min_in, lat_min_in, lon_max_in - lon_min_in, lat_max_in - lat_min_in],'LineWidth',2.5)


% draw letters for lake names
text(15.50,  -16.30, 'ERA-Interim global re-analysis', 'color', 'k', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Albert
text(22.50,  -11.30, 'CCLM 12km'                     , 'color', 'k', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Albert
text(29.70,  -4.30, 'CCLM 2.8km'                    , 'color', 'k', 'Fontweight', 'bold', 'FontSize', 11, 'FontName', 'Times New Roman') % Kivu


% save figure
export_fig DEM_C1_Domain -transparent


end



% --------------------------------------------------------------------
% DEM of Lake Victoria Basin
% --------------------------------------------------------------------


if flag_LVB == 1
    

% definegrid boundaries
lat_min_LVB = -5;  % = 31 - 5.625 - 5;
lat_max_LVB = 3;  % = 31 + 5.625 + 5;
lon_min_LVB = 27; % = -3.5 - 6.875 - 5;
lon_max_LVB = 36.5;  % = -3.6 + 6.875 + 5;


% call import and plot function
X = mf_readhgt(lat_min_LVB:lat_max_LVB, lon_min_LVB:lon_max_LVB, 'merge', 'plot');


% save figure
export_fig DEM_C1_LVB -transparent


end


toc
