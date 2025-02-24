% --------------------------------------------------------------------
% function to rescale Percatanges [0 100] so that they can be plotted  
% on logarithmic scale (90% 99% 99.9% 99.99%) 
% --------------------------------------------------------------------


function [PCT_plot] = mf_PCT_plotscale(PCT)

% input to the function are percentage values

% the function outputs rescaled percentage values ready for plotting



% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% initialize PR
PCT_plot = PCT;



% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------


% rescale PCT values for plotting on log scale
PCT_plot = 1 ./ (100 - PCT);




end

