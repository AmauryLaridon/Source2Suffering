

% load excel 
[unwpp_data_step1] = readtable('WPP2019_MORT_F16_1_LIFE_EXPECTANCY_BY_AGE_BOTH_SEXES_orig');


% select only country info/region and relevant columns
unwpp_data_step2_country = unwpp_data_step1(ismember(unwpp_data_step1.Type,'Country/Area'),[3 8 11]);
unwpp_data_step2_region  = unwpp_data_step1(ismember(unwpp_data_step1.Type,{'World', 'Income Group', 'Region', 'SDG region', 'Subregion'}),[3 8 11]);


% sort table in alphabetical order according to country/region name
unwpp_data_step3_country = sortrows(unwpp_data_step2_country, 1);
unwpp_data_step3_region  = sortrows(unwpp_data_step2_region , 1);


% get country/region names as a numeric vector
[unwpp_data_step4_country, ~, unwpp_countryindex] = unique(unwpp_data_step3_country(:,1));
[unwpp_data_step4_region , ~, unwpp_regionindex ] = unique(unwpp_data_step3_region( :,1));


% reshape table to fit world bank format
unwpp_data_step4_country{:,2:15} = cell2mat(accumarray(unwpp_countryindex, unwpp_data_step3_country.x5, [], @(x) {x(1:14).'}));
unwpp_data_step4_region{ :,2:15} = cell2mat(accumarray(unwpp_regionindex , unwpp_data_step3_region.x5 , [], @(x) {x(1:14).'}));


% save as excel table
writetable(unwpp_data_step4_country,'unwpp_data_step4_country.xlsx');
writetable(unwpp_data_step4_region,'unwpp_data_step4_region.xlsx');


% next step: manually transfer to world bank data format: 
% - take a world bank country file as example
% - sort data according to country
% - copy in unwpp_data_step4.xlsx
% - manually match the rows
% - sort by country code