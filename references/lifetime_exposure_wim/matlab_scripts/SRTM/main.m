

% --------------------------------------------------------------------
% DEM of selected domain
% --------------------------------------------------------------------


tic


% clean up
clc;
clear;
close all;


% flags
flag_DEM  = 1; % 0: do not plot SRTM DEM
               % 1: plot SRTM-DEM
flag_bigd = 1; % 0: do not plot big domain
               % 1: plot big domain
          
               

% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% add matlab scripts directory to path
addpath(genpath('C:\Users\u0079068\Documents\Research\matlab_scripts'));


% initialise domain grid parameters (from PEP file)
if flag_bigd == 0
    
    lat_min = -10.3750; % = -3.5 - 6.875;
    lat_max = 3.2750;   % = -3.6 + 6.875;
    lon_min = 25.3750;  % = 31 - 5.625;
    lon_max = 36.6250;  % = 31 + 5.625;
    
elseif flag_bigd == 1
    
    lat_min = -10.3750 - 5; % = -3.5 - 6.875 - 5;
    lat_max = 3.2750 + 5;   % = -3.6 + 6.875 + 5;
    lon_min = 25.3750 - 5;  % = 31 - 5.625 - 5;
    lon_max = 36.6250 + 5;  % = 31 + 5.625 + 5;
    
end



% --------------------------------------------------------------------
% DEM of selected domain
% --------------------------------------------------------------------


if flag_DEM == 1


% call import and plot function
X = mf_readhgt(lat_min:lat_max, lon_min:lon_max, 'merge', 'plot');



end


toc
