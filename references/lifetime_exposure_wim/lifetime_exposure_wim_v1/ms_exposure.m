
% --------------------------------------------------------------------
% subroutine to compute extreme event exposure across a person's lifetime
% note: preferably run "main"
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% convert Area Fraction Affected (AFA) to per-country 
% number of extremes affecting one individual across life span
% --------------------------------------------------------------------



% get exposure
if flags.exposure == 1   


    % prepare for loop
    exposure_peryear_percountry         = NaN(nruns, ncountries, nyears);
    exposure_peryear_perregion          = NaN(nruns, nregions  , nyears);
    
    landfrac_peryear_perregion          = NaN(nruns, nregions  , nyears);
    landfrac_peryear_perregion_15       = NaN(nruns, nregions  , nyears);
    landfrac_peryear_perregion_20       = NaN(nruns, nregions  , nyears);
    landfrac_peryear_perregion_NDC      = NaN(nruns, nregions  , nyears);
    landfrac_peryear_perregion_OS       = NaN(nruns, nregions  , nyears);
    landfrac_peryear_perregion_noOS     = NaN(nruns, nregions  , nyears);
    
    exposure_perrun_RCP                 = NaN(nruns, ncountries, nbirthyears);
    exposure_perrun_15                  = NaN(nruns, ncountries, nbirthyears);
    exposure_perrun_20                  = NaN(nruns, ncountries, nbirthyears);
    exposure_perrun_NDC                 = NaN(nruns, ncountries, nbirthyears);
    exposure_perrun_OS                  = NaN(nruns, ncountries, nbirthyears);
    exposure_perrun_noOS                = NaN(nruns, ncountries, nbirthyears);
    exposure_perrun_R26eval             = NaN(nruns, ncountries, nbirthyears);

    exposure_perregion_perrun_RCP       = NaN(nruns, nregions, nbirthyears);
    exposure_perregion_perrun_15        = NaN(nruns, nregions, nbirthyears);
    exposure_perregion_perrun_20        = NaN(nruns, nregions, nbirthyears);
    exposure_perregion_perrun_NDC       = NaN(nruns, nregions, nbirthyears);
    exposure_perregion_perrun_OS        = NaN(nruns, nregions, nbirthyears);  
    exposure_perregion_perrun_noOS      = NaN(nruns, nregions, nbirthyears);    
    exposure_perregion_perrun_R26eval   = NaN(nruns, nregions, nbirthyears);    
    

    if flags.embers == 1
    exposure_perregion_perrun_BE      = NaN(nruns, nregions, nbirthyears, nGMTsteps);
    end
    
    RCP2GMT_maxdiff_15                = NaN(nruns, 1);
    RCP2GMT_maxdiff_20                = NaN(nruns, 1);
    RCP2GMT_maxdiff_NDC               = NaN(nruns, 1);
    RCP2GMT_maxdiff_OS                = NaN(nruns, 1);
    RCP2GMT_maxdiff_noOS              = NaN(nruns, 1);
    RCP2GMT_maxdiff_R26eval           = NaN(nruns, 1);
    

    % loop over simulations
    for i=1:nruns
%     for i=250:261 % heatwaves
%     for i=250:nruns % heatwaves & TCs
    

        % print country name to screen
        disp(['run ' num2str(i) ' of ' num2str(nruns)])
        

        % Get ISIMIP GMT indices closest to GMT trajectories        
        [RCP2GMT_diff_15     , ind_RCP2GMT_15     ] = min(abs(bsxfun(@minus,isimip(i).GMT, GMT_15'  )));
        [RCP2GMT_diff_20     , ind_RCP2GMT_20     ] = min(abs(bsxfun(@minus,isimip(i).GMT, GMT_20'  )));
        [RCP2GMT_diff_NDC    , ind_RCP2GMT_NDC    ] = min(abs(bsxfun(@minus,isimip(i).GMT, GMT_NDC' )));
        [RCP2GMT_diff_OS     , ind_RCP2GMT_OS     ] = min(abs(bsxfun(@minus,isimip(i).GMT, GMT_OS'  )));
        [RCP2GMT_diff_noOS   , ind_RCP2GMT_noOS   ] = min(abs(bsxfun(@minus,isimip(i).GMT, GMT_noOS')));
        [RCP2GMT_diff_R26eval, ind_RCP2GMT_R26eval] = min(abs(bsxfun(@minus,isimip(i).GMT, isimip(1).GMT')));


        % Get maximum T difference between RCP and GMT trajectories (to remove rows later)
        RCP2GMT_maxdiff_15(i,1)      = nanmax(RCP2GMT_diff_15     );
        RCP2GMT_maxdiff_20(i,1)      = nanmax(RCP2GMT_diff_20     );
        RCP2GMT_maxdiff_NDC(i,1)     = nanmax(RCP2GMT_diff_NDC    );
        RCP2GMT_maxdiff_OS(i,1)      = nanmax(RCP2GMT_diff_OS     );
        RCP2GMT_maxdiff_noOS(i,1)    = nanmax(RCP2GMT_diff_noOS   );
        RCP2GMT_maxdiff_R26eval(i,1) = nanmax(RCP2GMT_diff_R26eval);


        % load AFA data of that run
        load(['ncfiles\workspaces\mw_isimip_AFA_' num2str(i)]);
        
        
        % --------------------------------------------------------------------
        % per country 
        % --------------------------------------------------------------------

        
        % get spatial average
        for j=1:ncountries 
            
            % corresponding picontrol - assume constant 1960 population density (this line takes about 16h by itself)
            [~, exposure_peryear_percountry_pic{i,j}] = mf_fieldmean(isimip_AFA_pic, population(:,:,1), countries.mask{j}); %#ok<SAGROW>
            
            % historical + RCP simulations
            [~, exposure_peryear_percountry(i,j,:)  ] = mf_fieldmean(isimip_AFA    , population       , countries.mask{j});
            
        end


        % call function to compute extreme event exposure per country and per lifetime
        exposure_perrun_RCP(i,:,:)     = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,:                  ) );
        exposure_perrun_15(i,:,:)      = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_15     ) );
        exposure_perrun_20(i,:,:)      = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_20     ) );
        exposure_perrun_NDC(i,:,:)     = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_NDC    ) );
        exposure_perrun_OS(i,:,:)      = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_OS     ) );
        exposure_perrun_noOS(i,:,:)    = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_noOS   ) );
        exposure_perrun_R26eval(i,:,:) = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_R26eval) );
        
        

        % --------------------------------------------------------------------
        % per region
        % --------------------------------------------------------------------


        % loop over regions
        for j=1:nregions 

            
            % get spatial average: historical + RCP simulations
            [~, landfrac_peryear_perregion(i,j,:)] = mf_fieldmean(isimip_AFA, grid_area, regions.mask{j});

            
            % get weighted spatial average
            exposure_perregion_perrun_RCP(  i,j,:)   = sum(exposure_perrun_RCP(    i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
            exposure_perregion_perrun_15(  i,j,:)    = sum(exposure_perrun_15(     i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
            exposure_perregion_perrun_20(  i,j,:)    = sum(exposure_perrun_20(     i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
            exposure_perregion_perrun_NDC( i,j,:)    = sum(exposure_perrun_NDC(    i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
            exposure_perregion_perrun_OS(  i,j,:)    = sum(exposure_perrun_OS(     i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
            exposure_perregion_perrun_noOS(i,j,:)    = sum(exposure_perrun_noOS(   i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
            exposure_perregion_perrun_R26eval(i,j,:) = sum(exposure_perrun_R26eval(i, regions.ind_member_countries{j, 1},:) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
       

        end

        
        % call function to compute land fraction exposed to extreme events per region and per year - for eulerian perspective figure
        landfrac_peryear_perregion_15(i,:,:)   = landfrac_peryear_perregion(i,:,ind_RCP2GMT_15  );
        landfrac_peryear_perregion_20(i,:,:)   = landfrac_peryear_perregion(i,:,ind_RCP2GMT_20  );
        landfrac_peryear_perregion_NDC(i,:,:)  = landfrac_peryear_perregion(i,:,ind_RCP2GMT_NDC );
        landfrac_peryear_perregion_OS(i,:,:)   = landfrac_peryear_perregion(i,:,ind_RCP2GMT_OS  );
        landfrac_peryear_perregion_noOS(i,:,:) = landfrac_peryear_perregion(i,:,ind_RCP2GMT_noOS);

        
        %  calculations for Burning Embers diagram
        if flags.embers == 1
            for l=1:nGMTsteps

                % Get ISIMIP GMT indices closest to GMT trajectories        
                [RCP2GMT_diff_BE(l,:), ind_RCP2GMT_BE(l,:) ] = min(abs(bsxfun(@minus,isimip(i).GMT, GMT_BE(:,l)' ))); %#ok<SAGROW>

                % Get maximum T difference between RCP and GMT trajectories (to remove rows later)
                RCP2GMT_maxdiff_BE(i,l) = nanmax(RCP2GMT_diff_BE(l,:));

                % call function to compute extreme event exposure per country and per lifetime - for burning embers diagram
                exposure_percountry_perrun_BE(i,:,:,l) = mf_exposure(isimip(i).years, birth_years, countries, exposure_peryear_percountry(i,:,ind_RCP2GMT_BE(l,:)) );
                
                
                % loop over regions
                for j=1:nregions 

                    % get weighted spatial average
                    exposure_perregion_perrun_BE(i,j,:,l) = sum(exposure_percountry_perrun_BE(i, regions.ind_member_countries{j, 1}, :, l) .* regions.cohort_weights{j}, 2, 'omitnan') ./ sum(regions.cohort_weights{j}, 2, 'omitnan');
        
                end
            end
        end



    end


    % save arrays as matlab workspace
    disp('saving mw_exposure')
    save('mw_exposure', 'RCP2GMT_maxdiff_15', 'RCP2GMT_maxdiff_20', 'RCP2GMT_maxdiff_NDC', 'RCP2GMT_maxdiff_OS', 'RCP2GMT_maxdiff_noOS', 'RCP2GMT_maxdiff_R26eval', ... 
                        'landfrac_peryear_perregion', 'landfrac_peryear_perregion_15', 'landfrac_peryear_perregion_20', 'landfrac_peryear_perregion_NDC'          , ...
                        'landfrac_peryear_perregion_OS', 'landfrac_peryear_perregion_noOS'                                                                        , ...
                        'exposure_perrun_RCP'   , 'exposure_perrun_15'              , 'exposure_perrun_20'              , 'exposure_perrun_NDC'                   , ...
                        'exposure_perrun_OS'    , 'exposure_perrun_noOS'            , 'exposure_perrun_R26eval'                                                   , ...
                        'exposure_perregion_perrun_RCP'   , 'exposure_perregion_perrun_15'  , 'exposure_perregion_perrun_20'     , 'exposure_perregion_perrun_NDC', ...
                        'exposure_perregion_perrun_OS'    , 'exposure_perregion_perrun_noOS', 'exposure_perregion_perrun_R26eval'                                 , ...
                        'exposure_peryear_percountry', 'exposure_peryear_percountry_pic'     , ...
                        'RCP2GMT_maxdiff_BE', 'exposure_percountry_perrun_BE', 'exposure_perregion_perrun_BE');

elseif flags.exposure == 0

    % load matlab workspace
    disp('loading mw_exposure')
    load('mw_exposure')
    
end



% --------------------------------------------------------------------
% Process picontrol data
% --------------------------------------------------------------------


% Define percentages for computing percentiles
percentages = [90 99.9 99.999]; % 1-in-10 to 1-in-100'000


% select only unique picontrol data sets
mod_gcm_ext          = cellfun(@(x,y,z) [x y z], {isimip.model}', {isimip.gcm}', {isimip.extreme}','un',0); % concatenate 'model', 'gcm' and 'extreme' into one string
[~, ind_pic_uniq, ~] = unique(mod_gcm_ext, 'stable');                                                       % get their unique values
isimip_pic           = isimip(ind_pic_uniq);                                                                 % subset isimip struct


% remove the alternative heatwave definitions from the selection prior to computing the mutli-model mean
if flags.coldwaves == 0
    % % ind_pic    = ind_pic_uniq(ind_pic_uniq & ~(ismember({isimip_pic.extreme}, 'heatwavedarea') & ~ismember({isimip_pic.model}, 'hwmid-humidex'))');      % the old one
    ind_pic    = ind_pic_uniq(ind_pic_uniq & ~(ismember({isimip_pic.extreme}, 'heatwavedarea') & ~ismember({isimip_pic.model}, 'hwmid99'))');                % the new one, proposed by Stefan
    % % ind_pic    = ind_pic_uniq(ind_pic_uniq & ~(ismember({isimip_pic.extreme}, 'heatwavedarea') & ~ismember({isimip_pic.model}, 'hwmid99-tasmax35'  ))'); % test sensitivity of BE plot
    % % ind_pic    = ind_pic_uniq(ind_pic_uniq & ~(ismember({isimip_pic.extreme}, 'heatwavedarea') & ~ismember({isimip_pic.model}, 'hwmid99-humidex40' ))'); % test sensitivity of BE plot
    % % ind_pic    = ind_pic_uniq(ind_pic_uniq & ~(ismember({isimip_pic.extreme}, 'heatwavedarea') & ~ismember({isimip_pic.model}, 'hwmid97p5-tasmax35'))'); % test sensitivity of BE plot
elseif flags.coldwaves == 1
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!! coldwave analysis !!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!');
    ind_pic    = ind_pic_uniq(ind_pic_uniq & ~(ismember({isimip_pic.extreme}, 'heatwavedarea') & ~ismember({isimip_pic.model}, 'cwmid99'  ))');            % run this one for cold wave plot!!!
end
isimip_pic = isimip(ind_pic);                                                                                                                                % subset isimip struct


% for small countries: replace NaN by vector of NaNs.
[nan_row, nan_col] = find(cell2mat(cellfun(@(x)any(isnan(x)),exposure_peryear_percountry_pic,'UniformOutput',false)));
for i=1:length(nan_row)
    exposure_peryear_percountry_pic{nan_row(i), nan_col(i)} = NaN(size(exposure_peryear_percountry_pic{nan_row(i), nan_col(i) - 1}));
end


% get picontrol exposure of 60-yr old - no duplicate pic runs - per country
[exposure_percountry_pic, exposure_percountry_pic_mean, ~, ~, ~] = mf_exposure_pic(isimip_pic, extremes, countries, percentages, [], exposure_peryear_percountry_pic(ind_pic, :));


% get picontrol exposure per age
if flags.exposure_pic == 1

    % quantify 'pure effect of increased life expectancy' from picontrol simulations - per country
    for ind_birthyear=1:nbirthyears
        
        % print country name to screen
        disp(['birth year ' num2str(ind_birthyear) ' of ' num2str(nbirthyears)])

        % compute picontrol exposure per country and per age
        [~ , exposure_percountry_pic_mean_perage(:,:,ind_birthyear) , ~, ~, ~ ] = mf_exposure_pic(isimip_pic, extremes, countries, percentages, ind_birthyear, exposure_peryear_percountry_pic(ind_pic, :) ); %#ok<SAGROW>
        
    end

    % save arrays as matlab workspace
    disp('saving mw_exposure_pic')
    if flags.coldwaves == 0
        save('mw_exposure_pic', 'exposure_percountry_pic_mean_perage');
    elseif flags.coldwaves == 1
        save('mw_exposure_pic_coldwaves', 'exposure_percountry_pic_mean_perage');
    end


elseif flags.exposure_pic == 0

    % load matlab workspace
    disp('loading mw_exposure_pic')
    if flags.coldwaves == 0
        load('mw_exposure_pic')
    elseif flags.coldwaves == 1
        load('mw_exposure_pic_coldwaves')
    end
    
end


% loop over regions
for ind_region=1:nregions 
    
    
    % get weighted spatial average of picontrol exposure !!! of 60-yr old !!! - no duplicate pic runs - per region
    exposure_perregion_pic_mean(:,ind_region) = sum(exposure_percountry_pic_mean(:, regions.ind_member_countries{ind_region, 1}) .* regions.cohort_weights{ind_region}(:,:,ages==age_ref), 2, 'omitnan') ./ sum(regions.cohort_weights{ind_region}(:,:,ages==age_ref), 2, 'omitnan');

    
    % get weighted spatial average of picontrol exposure per age to quantify 'pure effect of increased life expectancy' from picontrol simulations - per region
    exposure_perregion_pic_mean_perage(:,ind_region,:) = sum(exposure_percountry_pic_mean_perage(:, regions.ind_member_countries{ind_region, 1}, :) .* regions.cohort_weights{ind_region}(:,:,ages==age_ref), 2, 'omitnan') ./ sum(regions.cohort_weights{ind_region}(:,:,ages==age_ref), 2, 'omitnan');
    
    
    % get weighted spatial average of pooled exposure of 60-yr old from the bootstrapping - no duplicate pic runs - per region
    % loop over extremes
    for ind_extreme=1:nextremes+1
        
        % initialise empty array
        exposure_percountry_inregion_pic = [];
        
        % loop over countries in a region
        for ind_countryinregion=1:length(regions.member_countries{ind_region, 1})
            
            % get pooled exposure of 60-yr old from the bootstrapping and store in a matrix
            exposure_percountry_inregion_pic(:,ind_countryinregion) = exposure_percountry_pic{ind_extreme,ind_countryinregion}';
            
        end
        
        % compute weighted average of pooled exposure of 60-yr old
        exposure_perregion_pic{ind_extreme,ind_region} = permute(sum(exposure_percountry_inregion_pic .* regions.cohort_weights{ind_region}(:,:,ages==age_ref), 2, 'omitnan') ./ sum(regions.cohort_weights{ind_region}(:,:,ages==age_ref), 2, 'omitnan'), [2 1]);
        
    end
   
    
end



% --------------------------------------------------------------------
% compute averages across runs and sums across extremes 
% --------------------------------------------------------------------


% get indices corresponding to RCPs
ind_rcp26 = ismember({isimip.rcp}, 'rcp26'     )';
ind_rcp60 = ismember({isimip.rcp}, 'rcp60'     )';
ind_gfdl  = ismember({isimip.gcm}, 'gfdl-esm2m')';
ind_gmt   = contains({isimip.rcp}, 'rcp'       )';   % consider all indices

% remove the alternative heatwave definitions from the selection prior to computing the multi-model mean
if flags.coldwaves == 0
    % ind_gmt(ismember({isimip.extreme}, 'heatwavedarea') & ~ismember({isimip.model}, 'hwmid-humidex')) = 0;  % the old one
    ind_gmt(ismember({isimip.extreme}, 'heatwavedarea') & ~ismember({isimip.model}, 'hwmid99')) = 0;          % the new one, proposed by Stefan
elseif flags.coldwaves == 1
    ind_gmt(ismember({isimip.extreme}, 'heatwavedarea') & ~ismember({isimip.model}, 'cwmid99')) = 0;          % the new one, proposed by Stefan
end


% call function computing the multi-model mean (MMM) exposure and the Exposure Multiplication Factor (EMF)
[exposure_15            , exposure_mms_15            , EMF_15                                                                                                                            ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perrun_15               , [], RCP2GMT_maxdiff_15    , RCP2GMT_maxdiff_threshold);
[exposure_20            , exposure_mms_20            , EMF_20                                                                                                                            ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perrun_20               , [], RCP2GMT_maxdiff_20    , RCP2GMT_maxdiff_threshold);
[exposure_NDC           , exposure_mms_NDC           , EMF_NDC                                                                                                                           ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perrun_NDC              , [], RCP2GMT_maxdiff_NDC   , RCP2GMT_maxdiff_threshold);
[exposure_perregion_15  , exposure_perregion_mms_15  , EMF_perregion_15     , EMF_perregion_q25_15     , EMF_perregion_q75_15                                                            ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_15     , [], RCP2GMT_maxdiff_15    , RCP2GMT_maxdiff_threshold);
[exposure_perregion_20  , exposure_perregion_mms_20  , EMF_perregion_20     , EMF_perregion_q25_20     , EMF_perregion_q75_20                                                            ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_20     , [], RCP2GMT_maxdiff_20    , RCP2GMT_maxdiff_threshold);
[exposure_perregion_NDC , exposure_perregion_mms_NDC , EMF_perregion_NDC    , EMF_perregion_q25_NDC    , EMF_perregion_q75_NDC                                                           ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_NDC    , [], RCP2GMT_maxdiff_NDC   , RCP2GMT_maxdiff_threshold);
[exposure_perregion_OS  , exposure_perregion_mms_OS  , EMF_perregion_OS     , EMF_perregion_q25_OS     , EMF_perregion_q75_OS  , exposure_perregion_q25_OS  , exposure_perregion_q75_OS  ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_OS     , [], RCP2GMT_maxdiff_OS    , RCP2GMT_maxdiff_threshold);
[exposure_perregion_noOS, exposure_perregion_mms_noOS, EMF_perregion_noOS   , EMF_perregion_q25_noOS   , EMF_perregion_q75_noOS, exposure_perregion_q25_noOS, exposure_perregion_q75_noOS] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_noOS   , [], RCP2GMT_maxdiff_noOS  , RCP2GMT_maxdiff_threshold);



[exposure_perregion_R26    , ~                          , EMF_perregion_R26    , EMF_perregion_q25_R26    , EMF_perregion_q75_R26                                                           ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_rcp26 & ind_gfdl, ages, age_ref, exposure_perregion_perrun_RCP    , [], []                     , RCP2GMT_maxdiff_threshold);
[exposure_perregion_R26eval, ~                          , EMF_perregion_R26eval, EMF_perregion_q25_R26eval, EMF_perregion_q75_R26eval                                                       ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt   & ind_gfdl, ages, age_ref, exposure_perregion_perrun_R26eval, [], RCP2GMT_maxdiff_R26eval, RCP2GMT_maxdiff_threshold);



[~                      , ~                          , EMF_perregion_15_young2pic , ~                              , ~                                                                   ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_15  , exposure_perregion_pic_mean, RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);
[~                      , ~                          , EMF_perregion_20_young2pic , ~                              , ~                                                                   ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_20  , exposure_perregion_pic_mean, RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);
[~                      , ~                          , EMF_perregion_NDC_young2pic, EMF_perregion_q25_NDC_young2pic, EMF_perregion_q75_NDC_young2pic                                     ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perregion_perrun_NDC , exposure_perregion_pic_mean, RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);


% for ISIpedia article: country-level EMF youn2pic
[~                      , ~                          , EMF_15_young2pic , ~                              , ~                                                                             ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perrun_15            , exposure_percountry_pic_mean, RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);
[~                      , ~                          , EMF_20_young2pic , ~                              , ~                                                                             ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perrun_20            , exposure_percountry_pic_mean, RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);
[~                      , ~                          , EMF_NDC_young2pic, ~                              , ~                                                                             ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt  , ages, age_ref, exposure_perrun_NDC           , exposure_percountry_pic_mean, RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);



% loop over heatwave definitions - for heatwave sensitivity analysis requested by the reviewer
for ind_heatwave=1:length(model_names.heatwavedarea)

    
    % get indices corresponding to RCPs
    ind_gmt_heatwave(:, ind_heatwave) = contains({isimip.rcp}, 'rcp'  )';   % consider all indices
    
    
    % remove other heatwave definitions from the selection prior to computing the multi-model mean
    ind_gmt_heatwave(ismember({isimip.extreme}, 'heatwavedarea') & ~ismember({isimip.model}, lower(model_names.heatwavedarea{ind_heatwave})), ind_heatwave) = 0;          % the new one, proposed by Stefan
    
    % compute Exposure Multiplication Factor (EMF) per heatwave definition - per country
    [~, ~, EMF_heatwave_permodel_tmp_15 ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt_heatwave(:, ind_heatwave), ages, age_ref, exposure_perrun_15 , [], RCP2GMT_maxdiff_15  , RCP2GMT_maxdiff_threshold);
    [~, ~, EMF_heatwave_permodel_tmp_NDC] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt_heatwave(:, ind_heatwave), ages, age_ref, exposure_perrun_NDC, [], RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);
    
    
    % keep only the first row
    EMF_heatwave_permodel_15(      ind_heatwave,:,:) = EMF_heatwave_permodel_tmp_15(      5,:,:);
    EMF_heatwave_permodel_NDC(     ind_heatwave,:,:) = EMF_heatwave_permodel_tmp_NDC(     5,:,:);
    
    
    % compute Exposure Multiplication Factor (EMF) per heatwave definition - per region
    [exposure_perregion_heatwave_permodel_tmp_15,  ~, EMF_perregion_heatwave_permodel_tmp_15 , EMF_perregion_q25_heatwave_permodel_tmp_15  , EMF_perregion_q75_heatwave_permodel_tmp_15 ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt_heatwave(:, ind_heatwave), ages, age_ref, exposure_perregion_perrun_15 , [], RCP2GMT_maxdiff_15  , RCP2GMT_maxdiff_threshold);
    [exposure_perregion_heatwave_permodel_tmp_20,  ~, EMF_perregion_heatwave_permodel_tmp_20 , EMF_perregion_q25_heatwave_permodel_tmp_20  , EMF_perregion_q75_heatwave_permodel_tmp_20 ] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt_heatwave(:, ind_heatwave), ages, age_ref, exposure_perregion_perrun_20 , [], RCP2GMT_maxdiff_20  , RCP2GMT_maxdiff_threshold);
    [exposure_perregion_heatwave_permodel_tmp_NDC, ~, EMF_perregion_heatwave_permodel_tmp_NDC, EMF_perregion_q25_heatwave_permodel_tmp_NDC , EMF_perregion_q75_heatwave_permodel_tmp_NDC] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt_heatwave(:, ind_heatwave), ages, age_ref, exposure_perregion_perrun_NDC, [], RCP2GMT_maxdiff_NDC , RCP2GMT_maxdiff_threshold);
    
    
    % keep only the first row
    exposure_perregion_heatwave_permodel_15( ind_heatwave,:,:) = exposure_perregion_heatwave_permodel_tmp_15( 5,:,:);
    exposure_perregion_heatwave_permodel_20( ind_heatwave,:,:) = exposure_perregion_heatwave_permodel_tmp_20( 5,:,:);
    exposure_perregion_heatwave_permodel_NDC(ind_heatwave,:,:) = exposure_perregion_heatwave_permodel_tmp_NDC(5,:,:);
    EMF_perregion_heatwave_permodel_15(      ind_heatwave,:,:) = EMF_perregion_heatwave_permodel_tmp_15(      5,:,:);
    EMF_perregion_heatwave_permodel_20(      ind_heatwave,:,:) = EMF_perregion_heatwave_permodel_tmp_20(      5,:,:);
    EMF_perregion_heatwave_permodel_NDC(     ind_heatwave,:,:) = EMF_perregion_heatwave_permodel_tmp_NDC(     5,:,:);
    EMF_perregion_q25_heatwave_permodel_15(  ind_heatwave,:,:) = EMF_perregion_q25_heatwave_permodel_tmp_15(  5,:,:);
    EMF_perregion_q25_heatwave_permodel_20(  ind_heatwave,:,:) = EMF_perregion_q25_heatwave_permodel_tmp_20(  5,:,:);
    EMF_perregion_q25_heatwave_permodel_NDC( ind_heatwave,:,:) = EMF_perregion_q25_heatwave_permodel_tmp_NDC( 5,:,:);
    EMF_perregion_q75_heatwave_permodel_15(  ind_heatwave,:,:) = EMF_perregion_q75_heatwave_permodel_tmp_15(  5,:,:);
    EMF_perregion_q75_heatwave_permodel_20(  ind_heatwave,:,:) = EMF_perregion_q75_heatwave_permodel_tmp_20(  5,:,:);
    EMF_perregion_q75_heatwave_permodel_NDC( ind_heatwave,:,:) = EMF_perregion_q75_heatwave_permodel_tmp_NDC( 5,:,:);


end


%  calculations for Burning Embers diagram
% loop over GMT steps
for ind_BE=1:nGMTsteps


    % print GT level to screen
    disp(['GMT step ' num2str(ind_BE) ' of ' num2str(nGMTsteps)])


    % get multi-model mean exposure
    [exposure_perregion_BE(:,:,:,ind_BE), ~, EMF_perregion_young2pic_BE(:,:,:,ind_BE), ~, ~, ~, ~, EMF_perregion_young2pic_BE_harmean(:,:,:,ind_BE), EMF_perregion_young2pic_BE_geomexp(:,:,:,ind_BE)] = mf_exposure_mmm(extremes, {isimip.extreme}, ind_gmt, ages, age_ref, exposure_perregion_perrun_BE(:,:,:,ind_BE), exposure_perregion_pic_mean, RCP2GMT_maxdiff_BE(:,ind_BE), RCP2GMT_maxdiff_threshold); %#ok<SAGROW>

    
    % get multi-model mean EMF relative to 1960 birth cohort under the lowest scenario (constant present-day T, i.e. lower left corner of BE plot)
    EMF_perregion_young2ref_BE(:,:,:,ind_BE) = exposure_perregion_BE(:,:,:,ind_BE) ./ exposure_perregion_BE(:,:, ages==age_ref,1);


    % loop over extremes
    for i=1:nextremes+1

        % loop over regions
        for j=1:nregions
%             for j=nregions


            % compute inverse percentiles of BE exposure scenarios given pic exposure distribution
            PCT_perregion_young2pic_BE(i,j,:,ind_BE) = mf_invprctile(exposure_perregion_pic{i,j}', squeeze(exposure_perregion_BE(i,j,:,ind_BE)));


            % for verification - TBR
            % NOT YET WORKING COMPLETELY - IT DOES FOR MOST PERCENTILES BUT NOT FOR ALL !!!!!
            % TEST OTHER FUNCTION OPTIONS??
            VERIFY_INVPRCTILE(i,j,:,ind_BE) = prctile(exposure_perregion_pic{i,j}', squeeze(PCT_perregion_young2pic_BE(i,j,:,ind_BE)));


        end


    end


end




% % diagnostic plot for testing
% figure;
% for i=1:8; plot(isimip(i).years, isimip(i).GMT); hold on; end
% legend('gfdl rcp2.6', 'gfdl rcp6.0', 'hadgem rcp2.6', 'hadgem rcp6.0', 'ipsl rcp2.6', 'ipsl rcp6.0', 'miroc rcp2.6', 'miroc rcp6.0', 'location', 'northwest')

