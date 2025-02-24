% --------------------------------------------------------------------
% function to rescale Probability Ratios < 1 so that they can be plotted  
% on an intuïtive scale 
% --------------------------------------------------------------------


function [PR_plot, PR_plot_pos, PR_plot_neg] = mf_PR_plotscale(PR)

% input to the function are risk ratio values

% the function outputs rescaled risk ratio values ready for plotting



% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% initialize PR
PR_plot = PR;



% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------


% rescale PR values < 1
PR_plot(PR_plot < 1) = -1 ./ PR_plot(PR_plot < 1) + 2;


% store positive and negative indices separately
ind_pos              = find(PR_plot >= 1);
ind_neg              = find(PR_plot < 1);
PR_plot_pos          = NaN(size(PR_plot));
PR_plot_neg          = NaN(size(PR_plot));
PR_plot_pos(ind_pos) = PR_plot(ind_pos);
PR_plot_neg(ind_neg) = PR_plot(ind_neg);



end

