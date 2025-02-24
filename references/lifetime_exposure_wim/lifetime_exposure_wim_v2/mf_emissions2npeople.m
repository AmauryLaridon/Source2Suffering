

% --------------------------------------------------------------------
% function to compute additional number of people facing one additional
% climate extreme across their lifetime due to a certain amount of 
% CO2 emissions. 
% Also compute heat-related mortality associated with the emissions
% using the concept of 'mortality cost of carbon' (Bressler, 2021 NatComm)
% --------------------------------------------------------------------


function [nr_children_facing_extra_climate_extreme] = mf_emissions2npeople(CO2_emissions, TCRE, exposure_perregion_BE, birth_years, year_start, year_end, GMT_BE, valp_cohort_size_abs, rounding, ind_extreme)


% note: 
% - first calculation: nr of of people facing one additional climate extreme across their lifetime
% - second calculation: nr of heat-related deaths between today and 2100

% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% Manipulations: create bootstrap samples of life-accumulated exposure
% per picontrol simulation
% --------------------------------------------------------------------


% Compute change in GMT from emissions 
dGMT = TCRE * CO2_emissions;                            


% get number of birth years for loop over birth years
years_loop  = [year_end:-1:year_start]';
nbirthyears = length(years_loop);


% loop over birth years 2020 to 2010
for i=1:nbirthyears

    
    % extract the lifetime absolute number of heatwaves that a particular 
    % generation will experience under the burning embers pathways:
    valc_exposure_climate_extreme_newborns = squeeze(exposure_perregion_BE(ind_extreme,12, birth_years==years_loop(i),:)); % lifetime absolute climate extreme exposure for particular generation
    valc_GMT_2100                          = squeeze(GMT_BE(end,:));                                                       % 2100 GMT anomaly


    % fit a linear curve through the relation and get the slope
    valc_pf                             = polyfit(valc_GMT_2100,valc_exposure_climate_extreme_newborns,1);
    valc_slope_exposure_climate_extreme = valc_pf(1);


    % exctract the number of people in the cohort
    nr_newborns(i,1) = valp_cohort_size_abs(end-i+1,12);


    % Compute the average change in lifetime heatwave exposure for each birth cohort member
    nr_extra_climate_extremes_newborns(i,1) = valc_slope_exposure_climate_extreme * dGMT;
 

    % from the previous number, compute the number of children that are expected to face just one extra heatwave
    if     rounding == 0
        nr_children_facing_extra_climate_extreme(i,1) = floor(nr_newborns(i,1) * nr_extra_climate_extremes_newborns(i,1));                 % without rounding
    elseif rounding == 1 
        nr_children_facing_extra_climate_extreme(i,1) = floor(nr_newborns(i,1) * nr_extra_climate_extremes_newborns(i,1) ./ 1000) .* 1000; % with rounding to 1000s
    elseif rounding == 2 
        nr_children_facing_extra_climate_extreme(i,1) = floor(nr_newborns(i,1) * nr_extra_climate_extremes_newborns(i,1) ./ 100 ) .* 100;  % with rounding to 100s
    end
    nr_children_facing_extra_climate_extreme(i,1) =  max([nr_children_facing_extra_climate_extreme(i,1) 0]);                               % negative number people is an unphysical quantity, set cases when this occurs to zero
    
    
end


% get the sum across the years 2015-2020
nr_children_facing_extra_climate_extreme(end+1,1) = sum(nr_children_facing_extra_climate_extreme(1:end,1));



end

