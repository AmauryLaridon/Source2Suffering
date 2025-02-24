

% --------------------------------------------------------------------
% function to compute field maxima
% works for multiple masks now
% --------------------------------------------------------------------


function [var_secondmax] = mf_secondmax(var)



% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------

var_secondmax = nanmax(var(var<nanmax(var)));


end

