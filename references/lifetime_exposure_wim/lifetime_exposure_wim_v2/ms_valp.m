
% --------------------------------------------------------------------
% Compute values used in the paper
% note: preferably run "main"
% --------------------------------------------------------------------


clc;


% --------------------------------------------------------------------
% Abstract
% --------------------------------------------------------------------


% Range of EMFs for newborn under NDC
valp_EMF_range_newborn_NDC(1) = min(squeeze(EMF_perregion_NDC(1:nextremes, 12, ages == 0)));
valp_EMF_range_newborn_NDC(2) = max(squeeze(EMF_perregion_NDC(1:nextremes, 12, ages == 0)))



% --------------------------------------------------------------------
% Results
% --------------------------------------------------------------------


% land fraction exposed to heatwaves in 2020 and 2100
if flags.plot_fig1 == 1

valp_landfrac_15_2020  = landfrac_15_plot(years_SR15 == 2020)
valp_landfrac_20_2020  = landfrac_20_plot(years_SR15 == 2020)
valp_landfrac_NDC_2020 = landfrac_NDC_plot(years_SR15 == 2020)
    
valp_landfrac_15_2100  = landfrac_15_plot(years_SR15 == 2100)
valp_landfrac_NDC_2100 = landfrac_NDC_plot(years_SR15 == 2100)
    
end


% 60-yr old & newborn: nr heatwaves across lifetime
if flags.plot_fig1 == 1
    % nr of heatwaves
    valp_nr_heatwaves_60yrold     = round(mean(exposure_bars(1,:)))   % take the mean of 1.5°, 2° and NDC 
    valp_nr_heatwaves_0yrold      = round(exposure_bars(2,:))           

    % and their standard deviation
    valp_nr_heatwaves_60yrold_std = round(mean([exposure_perregion_mms_15( 5, 12, ages==age_ref) exposure_perregion_mms_20( 5, 12, ages==age_ref) exposure_perregion_mms_NDC( 5, 12, ages==age_ref)]))   % take the mean of 1.5°, 2° and NDC 
    valp_nr_heatwaves_0yrold_std  = round([exposure_perregion_mms_15( 5, 12, ages==age_young) exposure_perregion_mms_20( 5, 12, ages==age_young) exposure_perregion_mms_NDC( 5, 12, ages==age_young)])           
end


% 60-yr old & newborn: EMF (just read from the plot)


% 6-yr old at 3°C: 
valp_EMF_6yr_3deg_young2pic      = extremes_legend';
valp_EMF_6yr_3deg_young2pic(:,2) =  num2cell(round(squeeze(EMF_perregion_young2pic_BE(:, 12, ages == 6, GMT_steps == 3)), 0))


% 0-yr old at 3.5°C: 
valp_EMF_0yr_3p5deg_young2pic      = extremes_legend';
valp_EMF_0yr_3p5deg_young2pic(:,2) =  num2cell(round(squeeze(EMF_perregion_young2pic_BE(:, 12, ages == 0, GMT_steps == 3.5)), 0))


% aggregated exposure change of age cohorts below 20 under 1.5°
valp_EMF_allextremes_young2pic      = ages(ages<=20);
valp_EMF_allextremes_young2pic(:,2) = squeeze(EMF_perregion_young2pic_BE(nextremes + 1, 12, ages<=20, GMT_steps == 1.5));  % 1.5°C of warming
valp_EMF_allextremes_young2pic(:,3) = squeeze(EMF_perregion_young2pic_BE(nextremes + 1, 12, ages<=20, GMT_steps == 3))     % 3°C of warming


% EMF in Middle East and North Africa for age cohorts below 25: read from line plot


% regional differences in EMF for a newborn
valp_EMF_worldregions_newborns_NDC      = regions.name;
valp_EMF_worldregions_newborns_NDC(:,2) = num2cell(round(squeeze(EMF_perregion_NDC_young2pic(nextremes + 1, :, end)),1)')


% Change in burden per region (all extremes) - young2pic
valp_change_in_burden_perregion      = regions.name;
valp_change_in_burden_perregion_NDC  = exposure_perregion_NDC(nextremes+1, 1:12, ages == 0)' - exposure_perregion_pic_mean_perage(nextremes+1, 1:12, ages == age_ref)';
valp_change_in_burden_perregion_15   = exposure_perregion_15( nextremes+1, 1:12, ages == 0)' - exposure_perregion_pic_mean_perage(nextremes+1, 1:12, ages == age_ref)';
valp_change_in_burden_perregion(:,2) = num2cell(round((1 - valp_change_in_burden_perregion_15 ./ valp_change_in_burden_perregion_NDC) .* 100) .* -1)


% changes in relative cohort size
valp_relative_cohort_size_high_income_60 = round(regions.cohort_size_rel{3,1}   .* 100)
valp_relative_cohort_size_high_income_0  = round(regions.cohort_size_rel{3,end} .* 100)
valp_relative_cohort_size_low_income_60  = round(regions.cohort_size_rel{5,1}   .* 100)
valp_relative_cohort_size_low_income_0   = round(regions.cohort_size_rel{5,end} .* 100)


% number of children in EUCA and SSA born after 2015
% loop over ages
for ind_age=1:length(ages)
    
    % loop over regions
    for ind_region=1:nregions
        
        % get absolute cohort sizes of newborns for each region in 2020
        % transfer from thousands to actual numbers
        valp_cohort_size_abs(ind_age, ind_region)  = round(sum(regions.cohort_weights{ind_region}(:,:,ages==ages(ind_age))) .* 1000, 1);
        
        % get relative cohort sizes of newborns for each region in 2020
        valp_cohort_size_rel(ind_age, ind_region)  = round(sum(regions.cohort_weights{ind_region}(:,:,ages==ages(ind_age))) ./ sum(regions.cohort_weights{12}(:,:,ages==ages(ind_age))) .* 100);
        
    end
end
valp_cohort_size_abs_EUCA_0to5 = round(sum(valp_cohort_size_abs(end-5:end, 2 )), -6)
valp_cohort_size_abs_SSA_0to5  = round(sum(valp_cohort_size_abs(end-5:end, 10)), -6)


% EMF children under 5 - EUCA - young2pic
valp_EMF_EUCA_0to5_NDC_young2pic_heatwaves = squeeze(round(EMF_perregion_NDC_young2pic(5, 2, end-5:end), 1)) % not used
valp_EMF_EUCA_0to5_NDC_young2pic_all       = squeeze(round(EMF_perregion_NDC_young2pic(7, 2, end-5:end), 1)) % 3.8-4.0


% EMF children under 5 - SSA - young2pic
valp_EMF_SSA_0to5_NDC_young2pic_heatwaves  = squeeze(round(EMF_perregion_NDC_young2pic(5, 10, end-5:end), 1)) % 49-54
valp_EMF_SSA_0to5_NDC_young2pic_all        = squeeze(round(EMF_perregion_NDC_young2pic(7, 10, end-5:end), 1)) % 5.4-5.9


% change in globa life expectancy
valp_change_in_global_life_expectancy = [round(regions.life_expectancy_5{12}(1)) round(regions.life_expectancy_5{12}(end))]


% contribution of CC vs. LE change
if flags.plot_sfig11 == 1
disp('Percentage of CC signal compared to total change')

    for ind_region = 1:nregions
        
% %             % copied from ms_plotscript: prepare area plots
% %             exposure_perregion_percause_NDC = [repmat(exposure_perregion_pic_mean(nextremes + 1, ind_region), nbirthyears, 1) ...                                        % prepare for area plot
% %                                        flipud(squeeze(exposure_perregion_diff_cc( nextremes + 1, ind_region, :)))     ...
% %                                        flipud(squeeze(exposure_perregion_diff_dle(nextremes + 1, ind_region, :)))        ];
                            
            % copied from ms_plotscript: prepare area plots
            exposure_perregion_percause_NDC = [repmat(exposure_perregion_pic_mean_meanexp(1, ind_region), nbirthyears, 1) ...                                        % prepare for area plot
                                       flipud(squeeze(exposure_perregion_diff_cc(         1, ind_region, :)))     ...
                                       flipud(squeeze(exposure_perregion_diff_dle(        1, ind_region, :)))        ];

            % get percentage of CC signal compared to total change
            valp_CCvsLE_newborn(ind_region, 1) = round(exposure_perregion_percause_NDC(1,2) ./ (exposure_perregion_percause_NDC(1,2) + exposure_perregion_percause_NDC(1,3)) .* 100);
     
            % print to screen
            disp([regions.name{ind_region, 1}, ' ', num2str(valp_CCvsLE_newborn(ind_region, 1))])

     
    end
        
end


% Change in burden per extreme (global) - summary paragraph - young2pic
valp_change_in_burden_perextreme      = extremes_legend';
valp_change_in_burden_perextreme_NDC  = exposure_perregion_NDC(1:nextremes+1, 12, ages == 0) - exposure_perregion_pic_mean_perage(1:nextremes+1, 12, ages == age_ref);
valp_change_in_burden_perextreme_15   = exposure_perregion_15( 1:nextremes+1, 12, ages == 0) - exposure_perregion_pic_mean_perage(1:nextremes+1, 12, ages == age_ref);
valp_change_in_burden_perextreme(:,2) = num2cell(round((1 - valp_change_in_burden_perextreme_15 ./ valp_change_in_burden_perextreme_NDC) .* 100) .* -1)



% --------------------------------------------------------------------
% Methods
% --------------------------------------------------------------------


% average number of pic years per simulation:
nyears_pic_average = round(mean([isimip_pic.nyears_pic]))


% GMT anomalies in 2091-2100
valp_GMT_15_endofcentury  = round(mean(GMT_15( years_SR15 >= 2091 & years_SR15 <= 2100)),1)
valp_GMT_20_endofcentury  = round(mean(GMT_20( years_SR15 >= 2091 & years_SR15 <= 2100)),1)
valp_GMT_NDC_endofcentury = round(mean(GMT_NDC(years_SR15 >= 2091 & years_SR15 <= 2100)),1)



% --------------------------------------------------------------------
% Supplementary information
% --------------------------------------------------------------------


% ISIMIP simulation table:
% isimip_pic copied to excel and manually reworked
