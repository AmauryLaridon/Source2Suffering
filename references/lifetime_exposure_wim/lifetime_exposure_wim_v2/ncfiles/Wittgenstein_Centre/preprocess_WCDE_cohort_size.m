

% load excel 
[wcde_data_step1] = readtable('wcde_data_tmp');


% sort table in alphabetical order according to country/region name
wcde_data_step2 = sortrows(wcde_data_step1, 1);


% get country/region names as a numeric vector
[wcde_data_step3, ~, wcde_countryindex] = unique(wcde_data_step2(:,1));


% get number of columns
nagegroups = 21;
nyearblocks = 31;
ncols = nagegroups * nyearblocks;


% reshape table to fit world bank format
wcde_data_step3{:,2:ncols+1} = cell2mat(accumarray(wcde_countryindex, wcde_data_step2.Population, [], @(x) {x(1:ncols).'}));


% save as excel table
writetable(wcde_data_step3,'wcde_data_step3.xlsx');


% next step: manually transfer to world bank data format: 
% - take a world bank country excel file as example and remove the data
% - sort data according to country
% - copy in wcde_data_step4.xlsx
% - manually match the rows
% - sort by country code