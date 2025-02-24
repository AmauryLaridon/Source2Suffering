

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



% initialise circle parameters
lon_centre = 33; % corresponding to 33°E
lat_centre = -1; % corresponding to 1°S
circle_rad = 1.4;



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
rectangle('Position',[lon_min   , lat_min   , lon_max    - lon_min   , lat_max    - lat_min   ],'LineWidth',2.5); hold on;
% rectangle('Position',[lon_min_in, lat_min_in, lon_max_in - lon_min_in, lat_max_in - lat_min_in],'LineWidth',2.5)



% draw Lake Victoria circle
ang         = 0:0.01:2*pi; 
lon_circle  = circle_rad * cos(ang);
lat_circle  = circle_rad * sin(ang);
plot(lon_centre + lon_circle, lat_centre + lat_circle, 'color', [0.89 0.10 0.11], 'linewidth', 2);


% % draw transect line 29.5-36E; 1.5S
% line([29.5 36], [-1.50 -1.50], 'Color', 'k', 'LineWidth', 2);


% draw letters for lake names
text(33.0 , -3.00, 'a', 'Fontweight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman') % Victoria


% draw rectangle for the land mask ('surrounding land')
lat_min_Vict = -4.01;
lat_max_Vict = 2.01;
lon_min_Vict = 29.98;
lon_max_Vict = 35.95;
rectangle('Position',[lon_min_Vict, lat_min_Vict, lon_max_Vict - lon_min_Vict, lat_max_Vict - lat_min_Vict], 'EdgeColor', [0.89 0.10 0.11], 'LineWidth',2.5)


% save figure
export_fig DEM_CCLM_future -transparent

toc
