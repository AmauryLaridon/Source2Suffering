

% --------------------------------------------------------------------
% DEM of selected domain
% --------------------------------------------------------------------


tic


% clean up
clc;
clear;
close all;
          
               

% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% add matlab scripts directory to path
addpath(genpath('C:\Users\u0079068\Documents\Research\matlab_scripts'));



% --------------------------------------------------------------------
% DEM of selected domain
% --------------------------------------------------------------------


% definegrid boundaries
lat_min_EAfr = -10.3750 - 5; % = -3.5 - 6.875 - 5;
lat_max_EAfr = 3.2750 + 5;   % = -3.6 + 6.875 + 5;
lon_min_EAfr = 25.3750 - 5;  % = 31 - 5.625 - 5;
lon_max_EAfr = 36.6250 + 5;  % = 31 + 5.625 + 5;


% call import and plot function
X = mf_readhgt(lat_min_EAfr:lat_max_EAfr, lon_min_EAfr:lon_max_EAfr, 'merge', 'plot');


% draw inner and outer domain borders
dx         = 0.0625;
lon_min    = 31 - 5.625;
lon_max    = lon_min + dx * 180;
lat_min    = 1 - 11.375;
lat_max    = lat_min + dx * 220;
lon_min_in = lon_min + 10 * dx;
lon_max_in = lon_max - 10 * dx;
lat_min_in = lat_min + 10 * dx;
lat_max_in = lat_max - 10 * dx;
rectangle('Position',[lon_min   , lat_min   , lon_max    - lon_min   , lat_max    - lat_min   ],'LineWidth',2.5)
% rectangle('Position',[lon_min_in, lat_min_in, lon_max_in - lon_min_in, lat_max_in - lat_min_in],'LineWidth',2.5)


% draw transect line 26-36E; 0.55S
line([26 36], [-0.55 -0.55], 'Color', 'k', 'LineWidth', 2);


% draw letters for lake names
text(33.0 , -3.00, 'a', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Victoria
% text(30.0 , -5   , 'b', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Tanganyika
% text(30   ,  2.10, 'c', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Albert
% text(29.75, -2.00, 'd', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Kivu


% save figure
export_fig DEM_CCLM_present -transparent

toc
