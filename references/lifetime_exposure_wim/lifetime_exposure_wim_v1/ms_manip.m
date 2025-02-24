
% --------------------------------------------------------------------
% subroutine to perform manipulations on loaded variables
% note: preferably run "main"
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% manipulations: general
% --------------------------------------------------------------------


% get grid dimensions
[nlat, nlon] = size(lat_mod);


% get country names
countries.abbreviation = worldbank_country_meta(:,2);
countries.name         = worldbank_country_meta(:,1);
countries.region       = worldbank_country_meta(:,3);
countries.incomegroup  = worldbank_country_meta(:,4);


% get country median age in 2015 and 2020
countries.median_age_2015 = unwpp_medianage_country_data(:,1);
countries.median_age_2020 = unwpp_medianage_country_data(:,2);


% get country HDI
countries.hdi_2015 = unhdi_country_data(:,26);
countries.hdi_2018 = unhdi_country_data(:,29);


% get country HDI
countries.popunder5_2015 = un_popunder5_country_data(:,10);
countries.popunder5_2018 = un_popunder5_country_data(:,13);


% get region names
regions.abbreviation = worldbank_region_meta(:,2);
regions.name         = worldbank_region_meta(:,1);


% loop over regions to create 2-line version of region names (for pie charts)
for i=1:length(regions.name)
    if contains(regions.name{i}, ' & ')
        regions.name_2lines{i,1} = strrep(regions.name{i},' & ', ' \newline& ');    
    elseif contains(regions.name{i}, ' middle income')
        regions.name_2lines{i,1} = strrep(regions.name{i},' middle income', ' \newlinemiddle \newlineincome');    
    else 
        regions.name_2lines{i,1} = regions.name{i};
    end    
end


% loop over regions to create short version of region names (for pie charts)
for i=1:length(regions.name)
    if contains(regions.name{i}, ' middle income')
        regions.name_short{i,1} = strrep(regions.name{i},' middle income', ' \newlinemiddle');    
    elseif contains(regions.name{i}, ' income')
        regions.name_short{i,1} = strrep(regions.name{i},' income', '');    
    elseif contains(regions.name{i}, 'East Asia & Pacific')
        regions.name_short{i,1} = 'EASP';    
    elseif contains(regions.name{i}, 'Europe & Central Asia')
        regions.name_short{i,1} = 'EUCA';    
    elseif contains(regions.name{i}, 'Latin America & Caribbean')
        regions.name_short{i,1} = 'LAMC';    
    elseif contains(regions.name{i}, 'Middle East & North Africa')
        regions.name_short{i,1} = 'MENA';    
    elseif contains(regions.name{i}, 'North America')
        regions.name_short{i,1} = 'NAM';    
    elseif contains(regions.name{i}, 'South Asia')
        regions.name_short{i,1} = 'SAS';    
    elseif contains(regions.name{i}, 'Sub-Saharan Africa')
        regions.name_short{i,1} = 'SSA';    
    else 
        regions.name_short{i,1} = regions.name{i};
    end    
end



% get number of simulations
nruns = length(isimip);


% get number of birth years
nbirthyears = length(birth_years);


% get number of extreme impact categories
nextremes = length(extremes);


% get number of years considered in this study
nyears = length(year_start:year_end);



% --------------------------------------------------------------------
% manipulations: get country masks and interpolate life expectancies
% --------------------------------------------------------------------


% original data runs from 1960 to 2017 but we want estimates from 1960 to 2020
% add three columns of NaNs
worldbank_country_data(:,length(worldbank_years)+1:length(worldbank_years)+year_ref-worldbank_years(end)) = NaN;    
worldbank_region_data( :,length(worldbank_years)+1:length(worldbank_years)+year_ref-worldbank_years(end)) = NaN;    


% declare and populate interpolated World Bank/UN WPP data array; will be used for gap filling
worldbank_country_data_interp = worldbank_country_data;
worldbank_region_data_interp  = worldbank_region_data;
unwpp_country_data_interp     = NaN(size(worldbank_country_data_interp));
unwpp_region_data_interp      = NaN(size(worldbank_region_data_interp));


% initialise index denoting too small countries or empty regions
ind_small_countries = [];
ind_empty_regions   = [];


% get country masks
if flags.masks == 1   

    % loop over countries
    for i=1:length(countries.name)


        % print country name to screen
        disp(countries.name{i})


        % store birth_year data
        countries.birth_years{i,1} = birth_years;


        % extract life expectancy at birth data from World Bank file and fill in missing data - not used in final analysis
        ind_nan                                  = find(isnan(worldbank_country_data(i,:)));
        ind_data                                 = find(~isnan(worldbank_country_data(i,:)));
        worldbank_country_data_interp(i,ind_nan) = interp1(countries.birth_years{i,1}(ind_data), worldbank_country_data(i,ind_data), countries.birth_years{i,1}(ind_nan), 'linear', 'extrap');
        countries.life_expectancy_0{i,1}         = worldbank_country_data_interp(i,:);


        % extract life expectancy at age 5 data from UN WPP file and
        % linearly interpolate from 5-year WPP blocks to pre-defined birth
        % year (extrapolate from 2013 to 2020, note that UN WPP has no NaNs)
        ind_data                         = find(~isnan(unwpp_country_data(i,:)));
        unwpp_country_data_interp(i,:)   = interp1(unwpp_years, unwpp_country_data(i,ind_data), countries.birth_years{i,1}, 'linear', 'extrap');
        countries.life_expectancy_5{i,1} = unwpp_country_data_interp(i,:) + 5 + 6; % add 5 to transfer from 'expected years left to live for 5-year old' to 'life expectancy of birth cohort excluding infant mortality'; add 6 to transfer from 'period life expectancy' to 'cohort life expectancy' as requested by Reviewer 1 (suggested by Marina based on Goldstein paper) (note that calendar years were already corrected during loading)


        % extract population size per age cohort data from WCDE file and
        % linearly interpolate from 5-year WCDE blocks to pre-defined birth year
        wcde_country_data_reshape   = reshape(wcde_country_data(i,:), length(wcde_years), []);                % reshape vector to matrix
        wcde_country_data_reshape   = [wcde_country_data_reshape(:,1), wcde_country_data_reshape];            % assume constant cohort sizes beyond array bounds (age 0 to 2)
        wcde_country_data_reshape   = [wcde_country_data_reshape; wcde_country_data_reshape(end,:)];          % assume constant cohort sizes beyond array bounds (years 2100 to 2113)       
        [Xorig, Yorig]              = meshgrid([min(ages); wcde_ages], [wcde_years; max(years_SR15)]);        % prepare for 2D interpolation
        [Xnew, Ynew]                = meshgrid(ages, years_SR15);                                             % prepare for 2D interpolation
        wcde_country_data_interp    = interp2(Xorig, Yorig, wcde_country_data_reshape, Xnew, Ynew, 'linear'); % 2D interpolation
        countries.cohort_size{i,1}  = wcde_country_data_interp ./ 5;                                          % divide by 5 to go from 5-year cohort size to 1-year cohort size


        % extract data from country borders file
        ind_borders            = find(strcmp(country_borders.dbf.NAME_LONG(:,1), countries.name{i}));
        countries.mask{i,1}    = inpolygon(lon_mod, lat_mod, country_borders.ncst{ind_borders}(:,1), country_borders.ncst{ind_borders}(:,2) );
        countries.borders{i,1} = country_borders.ncst{ind_borders};


        % extract total population
        countries.population{i,1} = squeeze(nansum( population .* repmat(countries.mask{i,1}, [1 1 nyears]), [1 2]));

        
        % check whether country has non-zero mask
        if isempty(find(countries.mask{i}==1, 1))
            ind_small_countries = [ind_small_countries i]; %#ok<AGROW>
        end


    end

    
    % Keep only those countries with non-zero masks
    countries.abbreviation(ind_small_countries)      = [];
    countries.name(ind_small_countries)              = [];
    countries.region(ind_small_countries)            = [];
    countries.incomegroup(ind_small_countries)       = [];
    countries.birth_years(ind_small_countries)       = [];
    countries.life_expectancy_0(ind_small_countries) = [];
    countries.life_expectancy_5(ind_small_countries) = [];
    countries.cohort_size(ind_small_countries)       = [];
    countries.population(ind_small_countries)        = [];
    countries.mask(ind_small_countries)              = [];
    countries.borders(ind_small_countries)           = [];
    countries.median_age_2015(ind_small_countries)   = [];
    countries.median_age_2020(ind_small_countries)   = [];
    countries.hdi_2015(ind_small_countries)          = [];
    countries.hdi_2018(ind_small_countries)          = [];
    countries.popunder5_2015(ind_small_countries)    = [];
    countries.popunder5_2018(ind_small_countries)    = [];

 
    % loop over regions
    for i=1:length(regions.name)


        % print region name to screen
        disp(regions.name{i})


        % store birth_year data
        regions.birth_years{i,1} = birth_years;


        % extract data from World Bank file and fill in missing data
        ind_nan                                 = find(isnan(worldbank_region_data(i,:)));
        ind_data                                = find(~isnan(worldbank_region_data(i,:)));
        worldbank_region_data_interp(i,ind_nan) = interp1(regions.birth_years{i,1}(ind_data), worldbank_region_data(i,ind_data), regions.birth_years{i,1}(ind_nan), 'linear', 'extrap');
        regions.life_expectancy_0{i,1}          = worldbank_region_data_interp(i,:);
        
        
        % extract life expectancy at age 5 data from UN WPP file and
        % linearly interpolate from 5-year WPP blocks to pre-defined birth
        % year (extrapolate from 2013 to 2020, note that UN WPP has no NaNs)
        ind_data                         = find(~isnan(unwpp_region_data(i,:)));
        unwpp_region_data_interp(i,:)    = interp1(unwpp_years, unwpp_region_data(i,ind_data), regions.birth_years{i,1}, 'linear', 'extrap');
        regions.life_expectancy_5{i,1}   = unwpp_region_data_interp(i,:) + 5 + 6; % add 5 to transfer from 'expected years left to live for 5-year old' to 'life expectancy of birth cohort excluding infant mortality'; add 6 to transfer from 'period life expectancy' to 'cohort life expectancy' as requested by Reviewer 1 (suggested by Marina based on Goldstein paper) (note that calendar years were already corrected during loading)

        
        % identify all countries belonging to the region
        ind_member_countries              = ismember(countries.region(:,1), regions.name{i}) | ismember(countries.incomegroup(:,1), regions.name{i});
        if strcmp(regions.name{i,1}, 'World')
            ind_member_countries = true(length(countries.name),1);
        end
        regions.ind_member_countries{i,1} = ind_member_countries;
        regions.member_countries{i,1}     = countries.name(ind_member_countries);

        
        % get total population in the region per cohort in 2020: 
        tmp1 = [];
        tmp2 = [];
        tmp1 = countries.cohort_size(ind_member_countries);
        for j=1:length(tmp1)
            tmp2(j,:) = tmp1{j,1}(years_SR15==year_ref, :);
        end
        regions.cohort_weights{i,1} = permute(tmp2, [3 1 2]);

        
        % merge country masks to generate the region mask
        tmp               = countries.mask(ind_member_countries);
        regions.mask{i,1} = logical(sum(cat(3, tmp{:}),3));

        
        % sum all masks to obtain the World mask
        if strcmp(regions.name{i,1}, 'World')
            regions.mask{i,1} = logical(sum(cat(3,countries.mask{:}),3));
        end
            

        % check whether region has non-zero mask
        if isempty(find(regions.mask{i}==1, 1))
            ind_empty_regions = [ind_empty_regions i]; %#ok<AGROW>
        end


    end
    
    
    % Keep only those regions with non-zero masks
    regions.abbreviation(ind_empty_regions)         = [];
    regions.name(ind_empty_regions)                 = [];
    regions.birth_years(ind_empty_regions)          = [];
    regions.life_expectancy_0(ind_empty_regions)    = [];
    regions.life_expectancy_5(ind_empty_regions)    = [];
    regions.ind_member_countries(ind_empty_regions) = [];
    regions.member_countries(ind_empty_regions)     = [];
    regions.cohort_weights(ind_empty_regions)       = [];
    regions.mask(ind_empty_regions)                 = [];
    
    
    % save struct as matlab workspace
    disp('saving mw_countries')
    save('mw_countries','countries', 'regions');


elseif flags.masks == 0

    % load matlab workspace
    disp('loading mw_countries')
    load('mw_countries')
    
end


% get number of countries and regions
ncountries = length(countries.name);
nregions   = length(regions.name);


% get total population of country in 2020
for i=1:ncountries
    population_percountry_2020(i,1) = countries.population{i,1}(years_pop == 2020); %#ok<SAGROW>
end
