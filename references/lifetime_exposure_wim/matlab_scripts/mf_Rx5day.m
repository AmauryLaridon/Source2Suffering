

% --------------------------------------------------------------------
% function to load CORDEX model data
% input: daily precipitation data and years in which they fall
% --------------------------------------------------------------------


% this function was inspired by this example
% https://meteo.unican.es/trac/MLToolbox/wiki/Observations/Extremes
% original function available here:
% https://meteo.unican.es/trac/MLToolbox/browser/MLToolbox/branches/MLToolbox_R2013/MeteoLab/Indicators/extremesIndicator.m


% definition of Rx5day from http://etccdi.pacificclimate.org/list_27_indices.shtml:
% Rx5day, ANNUAL maximum consecutive 5-day precipitation:
% Let RRkj be the precipitation amount for the 5-day interval ending k, period j. Then maximum 5-day values for period j are:
% Rx5dayj = max (RRkj)


function [Rx5day] = mf_Rx5day(TOT_PREC, year)


% get number of years
years  = (nanmin(year):nanmax(year))';
nyears = length(years);


% prepare for loop
Rx5day = NaN(nyears,size(TOT_PREC,2));


% loop over years
for i=1:nyears
    
    % get all observations in that year
	ind   = find(year == years(i))';
    
    % prepare for loop
	R5day = NaN(length(ind)-4,size(TOT_PREC,2));
    
    % loop over days in particular year
    for j = 1:length(ind)-4
         R5day(j,:) = nansum(TOT_PREC(ind(j):ind(j)+4,:)); % 5-day accumulated precipitation
    end
    
    % Annual maximum consecutive 5-day precipitation
    Rx5day(i,:) = nanmax(R5day);
    
end

end

