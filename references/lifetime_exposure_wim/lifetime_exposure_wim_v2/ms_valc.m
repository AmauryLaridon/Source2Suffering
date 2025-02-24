
% --------------------------------------------------------------------
% Generate numbers and Excel files for 
% Save The Children Report 
% Carbon Brief
% IPCC AR5
% The Guardian
% note: preferably run "main"
% --------------------------------------------------------------------


% flags for valSTC calculations
flags.valSTC_00   = 0; % 0: do not store numbers used in the report
                       % 1: store numbers used in the report
flags.valSTC_01   = 0; % 0: do not store birth cohort size
                       % 1: store birth cohort size
flags.valSTC_02   = 0; % 0: do not store EMF young2pic excel sheets
                       % 1: store EMF young2pic excel sheets
flags.valSTC_03   = 0; % 0: do not store EMF young2ref excel sheets
                       % 1: store EMF young2ref excel sheets
flags.valSTC_04   = 0; % 0: do not store EMF young2ref excel sheets per country
                       % 1: store EMF young2ref excel sheets per country
flags.valSTC_05   = 0; % 0: do not store EMF lowincome2highincome excel sheets
                       % 1: store EMF lowincome2highincome excel sheets
flags.valSTC_06   = 0; % 0: do not store EMF young2pic global allextremes excel sheets
                       % 1: store EMF young2pic global allextremes excel sheets
flags.valSTC_07   = 0; % 0: do not store EMF young2pic excel sheets per country
                       % 1: store EMF young2pic excel sheets per country
flags.valNHRI_01  = 0; % 0: do not store EMF young2pic excel sheets for Portugal, Duarte Age groups
                       % 1: store EMF young2pic excel sheets for Portugal, Duarte Age groups
flags.valNHRI_02  = 0; % 0: do not store reduced burden young2pic excel sheets for Portugal, Duarte Age groups
                       % 1: store EMF reduced burden excel sheets for Portugal, Duarte Age groups
flags.valVerus_01 = 0; % 0: do not store EMF young2pic excel sheets for Zwitserland, Duarte Age groups
                       % 1: store EMF young2pic excel sheets for Zwitserland, Duarte Age groups
flags.valVerus_02 = 0; % 0: do not store reduced burden young2pic excel sheets for Zwitserland, Duarte Age groups
                       % 1: store EMF reduced burden excel sheets for Zwitserland, Duarte Age groups
flags.valCB_01    = 0; % 0: do not store EMF young2pic excel sheets
                       % 1: store EMF young2pic excel sheets
flags.valCB_02    = 0; % 0: do not store EMF young2pic excel sheets
                       % 1: store EMF young2pic excel sheets
flags.valCB_03    = 0; % 0: do not store EMF young2pic excel sheets
                       % 1: store EMF young2pic excel sheets
flags.valAR6      = 0; % 0: do not store numbers for Europe
                       % 1: store numbers for Europe
flags.valBBC      = 0; % 0: do not store numbers for The BBC
                       % 1: store numbers for The BBC
flags.valTG       = 0; % 0: do not store numbers for The Guardian
                       % 1: store numbers for The Guardian
flags.valDPA      = 0; % 0: do not store numbers for DPA (German press agency)
                       % 1: store numbers for DPA (German press agency)
flags.valWP       = 0; % 0: do not store numbers for The Washington Post
                       % 1: store numbers for The Washington Post
flags.valISIpedia = 0; % 0: do not store numbers for ISIpedia
                       % 1: store numbers for ISIpedia
flags.valMCF      = 0; % 0: do not store numbers for MyClimateFuture
                       % 1: store numbers for MyClimateFuture
flags.valsciam    = 0; % 0: do not store numbers for scientific Amrican
                       % 1: store numbers for scientific Amrican
flags.valBloom    = 0; % 0: do not store numbers for Bloom
                       % 1: store numbers for Bloom
flags.valNorLic   = 0; % 0: do not store numbers for Norwegian lawsuit
                       % 1: store numbers for Norwegian lawsuit
flags.valNorStop  = 0; % 0: do not store numbers for Norwegian lawsuit - written report April 2024 for written proceedings
                       % 1: store numbers for Norwegian lawsuit
flags.valNorAppea = 0; % 0: do not store numbers for Norwegian lawsuit - written report June 2024 for appeal
                       % 1: store numbers for Norwegian lawsuit
flags.valNorECtHR = 0; % 0: do not store numbers for Norwegian lawsuit - written report June 2024 for appeal
                       % 1: store numbers for Norwegian lawsuit
flags.valUKrep1   = 0; % 0: do not store numbers for UK lawsuit - first written report october 2024
                       % 1: store numbers for UK lawsuit
flags.valUKrep2   = 0; % 0: do not store numbers for UK lawsuit - second written report october 2024
                       % 1: store numbers for UK lawsuit
flags.valRomania  = 1; % 0: do not store numbers for Romania communication
                       % 1: store numbers for Romania communication



% clean up
clc                       



% --------------------------------------------------------------------
% numbers used in the report
% --------------------------------------------------------------------
                       
                       
if flags.valSTC_00 == 1


% Change in exposure per extreme (global) 
valc_EMF_young2ref_perextreme      = extremes_legend';
valc_EMF_young2ref_perextreme(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(1:nextremes+1, 12, ages == 0)), 1))
                       

% Change in exposure per region (wildfires) 
valc_EMF_young2ref_perregion_NDC_wildfires      = regions.name;
valc_EMF_young2ref_perregion_NDC_wildfires(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(1, :, ages == 0)), 1))

                       
% Change in exposure per region (crop failure) 
valc_EMF_young2ref_perregion_NDC_cropfailures      = regions.name;
valc_EMF_young2ref_perregion_NDC_cropfailures(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(2, :, ages == 0)), 1))


% Change in exposure per region (droughts) 
valc_EMF_young2ref_perregion_NDC_droughts      = regions.name;
valc_EMF_young2ref_perregion_NDC_droughts(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(3, :, ages == 0)), 1))


% Change in exposure per region (river floods) 
valc_EMF_young2ref_perregion_NDC_riverfloods      = regions.name;
valc_EMF_young2ref_perregion_NDC_riverfloods(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(4, :, ages == 0)), 1))


% Change in exposure per region (heatwvaes) 
valc_EMF_young2ref_perregion_NDC_heatwaves      = regions.name;
valc_EMF_young2ref_perregion_NDC_heatwaves(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(5, :, ages == 0)), 1))


% Change in exposure per region (tropical cyclones) 
valc_EMF_young2ref_perregion_NDC_tropicalcyclones      = regions.name;
valc_EMF_young2ref_perregion_NDC_tropicalcyclones(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(6, :, ages == 0)), 1))


% change in INTERGENERATIONAL BURDEN per extreme (global) - summary paragraph - young2ref
valc_change_in_burden_young2ref_perextreme      = extremes_legend';
valc_change_in_burden_young2ref_perextreme_NDC  = exposure_perregion_NDC(1:nextremes+1, 12, ages == 0) - exposure_perregion_NDC(1:nextremes+1, 12, ages == age_ref);
valc_change_in_burden_young2ref_perextreme_15   = exposure_perregion_15( 1:nextremes+1, 12, ages == 0) - exposure_perregion_15( 1:nextremes+1, 12, ages == age_ref);
valc_change_in_burden_young2ref_perextreme(:,2) = num2cell(round((1 - valc_change_in_burden_young2ref_perextreme_15 ./ valc_change_in_burden_young2ref_perextreme_NDC) .* 100) .* -1)


% REDUCTION in LIFETIME EXPOSURE OF NEWBORNS per extreme (global) when following 1.5°C instead of NDC - STC EXTRA REQUEST 16/09/2021 - young2ref
valc_change_in_burden_newborns_NDC215_perextreme      = extremes_legend';
valc_change_in_burden_newborns_NDC215_perextreme(:,2) = num2cell(round(1 ./ (exposure_perregion_15(1:nextremes+1, 12, ages == 0) ./ exposure_perregion_NDC( 1:nextremes+1, 12, ages == 0)), 2 )); % EMF
valc_change_in_burden_newborns_NDC215_perextreme(:,3) = num2cell(round((exposure_perregion_NDC(1:nextremes+1, 12, ages == 0) - exposure_perregion_15(1:nextremes+1, 12, ages == 0)) ./ exposure_perregion_NDC( 1:nextremes+1, 12, ages == 0) .* 100 .* -1))  % percentage reduction


% INCREASE in LIFETIME EXPOSURE OF NEWBORNS per extreme (global) when following NDC instead of 1.5°C - STC EXTRA REQUEST 16/09/2021 - young2ref
valc_change_in_burden_newborns_152NDC_perextreme      = extremes_legend';
valc_change_in_burden_newborns_152NDC_perextreme(:,2) = num2cell(round((exposure_perregion_NDC(1:nextremes+1, 12, ages == 0)                                                         ./ exposure_perregion_15( 1:nextremes+1, 12, ages == 0))        , 2 )); % EMF
valc_change_in_burden_newborns_152NDC_perextreme(:,3) = num2cell(round((exposure_perregion_15( 1:nextremes+1, 12, ages == 0) - exposure_perregion_NDC(1:nextremes+1, 12, ages == 0)) ./ exposure_perregion_15( 1:nextremes+1, 12, ages == 0) .* 100 .* -1))        % percentage reduction


end



% --------------------------------------------------------------------
% 01 birth cohort size
% --------------------------------------------------------------------


if flags.valSTC_01 == 1


% loop over ages
for ind_age=1:length(ages)

        f
    % store age
    matrix_out_01_abs{ind_age+1, 1} = ages(ind_age);
    matrix_out_01_rel{ind_age+1, 1} = ages(ind_age);
    
    
    % loop over regions
    for ind_region=1:nregions

        
        % get absolute cohort sizes of newborns for each region in 2020
        % transfer from thousands to actual numbers
        matrix_out_01_abs{ind_age+1, ind_region+1}  = round(sum(regions.cohort_weights{ind_region}(:,:,ages==ages(ind_age))) .* 1000, -6);

        
        % get relative cohort sizes of newborns for each region in 2020
        matrix_out_01_rel{ind_age+1, ind_region+1}  = round(sum(regions.cohort_weights{ind_region}(:,:,ages==ages(ind_age))) ./ sum(regions.cohort_weights{12}(:,:,ages==ages(ind_age))) .* 100);

        
        % store region name
        matrix_out_01_abs{1, ind_region+1} = regions.name{ind_region};
        matrix_out_01_rel{1, ind_region+1} = regions.name{ind_region};
        
        
    end

    
end
        

% save as excell file
filename = 'follow-up_STC\data\01_cohort_size.xlsx';
writecell(matrix_out_01_abs, filename, 'Sheet', 'absolute numbers');



end



% --------------------------------------------------------------------
% 02 EMF young2pic excel sheets
% --------------------------------------------------------------------


if flags.valSTC_02 == 1


% loop over regions
for ind_region = 1:12
% for ind_region = 12
  
    
    % loop over extreme event categories
    for ind_extreme = 1:nextremes+1
%     for ind_extreme = nextremes+1

        
        % get associated Exposure Multiplication Factors (EMF)
        EMF_region_plot_BE = round(squeeze(EMF_perregion_young2pic_BE(ind_extreme, ind_region, :, :))', 1);

        
        % prepare output matrix
        matrix_out_02                     = NaN(nGMTsteps+1, nbirthyears+1);
        matrix_out_02(2:nGMTsteps+1, 1)   = fliplr(round(GMT_BE(end,:), 1));
        matrix_out_02(1, 2:nbirthyears+1) = ages;
        matrix_out_02(2:end, 2:end)       = flipud(EMF_region_plot_BE); 
        
        
        % save as excell file
        filename = ['follow-up_STC\data\02_EMF_young2pic_' regions.name{ind_region} '.xlsx'];
        writematrix(matrix_out_02, filename, 'Sheet', extremes_legend{ind_extreme});
 
        
    end

         
    % change cell color of first row and column to yellow
    excelObject = actxserver('Excel.Application');            % Launch an Excel server using ActiveX (Windows ONLY).
    fullFileName = fullfile(pwd, filename);                   % Create the filename of the existing workbook.
    excelWorkbook = excelObject.workbooks.Open(fullFileName); % Open the workbook from disk.
    excelObject.Visible = true;                               % Excel is invisible so far.  Make it visible.
    Excel_utils.FormatCellColor(excelWorkbook, 'A1..A29', 6); % WT change cell color of first column to yellow
    Excel_utils.FormatCellColor(excelWorkbook, 'A1..BJ1', 6); % WT change cell color of first row to yellow
    excelWorkbook.Save;                                       % Save the current state of the workbook.
    excelWorkbook.Close;                                      % Close the workbook.  Excel will stay open but be hidden.
    excelObject.Quit;                                         % Shut down the Excel server instance.
    clear('excelObject', 'excelWorkbook', 'yourFormula');     % Clear the excel object variable from MATLAB's memory.
    fprintf('Done interacting with Excel.\n');                % The clear finally shuts down the server and it no longer appears in Task Manager.
    

end



end



% --------------------------------------------------------------------
% 03 EMF young2ref excel sheets
% --------------------------------------------------------------------


if flags.valSTC_03 == 1   


% loop over regions
for ind_region = 1:12
% for ind_region = 12
  
    
    % loop over extreme event categories
    for ind_extreme = 1:nextremes+1
%     for ind_extreme = nextremes+1

        
        % get associated Exposure Multiplication Factors (EMF)
        EMF_region_plot_BE = round(squeeze(EMF_perregion_young2ref_BE(ind_extreme, ind_region, :, :))', 1);

        
        % prepare output matrix
        matrix_out_03                     = NaN(nGMTsteps+1, nbirthyears+1);
        matrix_out_03(2:nGMTsteps+1, 1)   = fliplr(round(GMT_BE(end,:), 1));
        matrix_out_03(1, 2:nbirthyears+1) = ages;
        matrix_out_03(2:end, 2:end)       = flipud(EMF_region_plot_BE); 
        
        
        % save as excell file
        filename = ['follow-up_STC\data\03_EMF_young2ref_' regions.name{ind_region} '.xlsx'];
        writematrix(matrix_out_03, filename, 'Sheet', extremes_legend{ind_extreme});
 
        
    end

         
    % change cell color of first row and column to yellow
    excelObject = actxserver('Excel.Application');            % Launch an Excel server using ActiveX (Windows ONLY).
    fullFileName = fullfile(pwd, filename);                   % Create the filename of the existing workbook.
    excelWorkbook = excelObject.workbooks.Open(fullFileName); % Open the workbook from disk.
    excelObject.Visible = true;                               % Excel is invisible so far.  Make it visible.
    Excel_utils.FormatCellColor(excelWorkbook, 'A1..A29', 6); % WT change cell color of first column to yellow
    Excel_utils.FormatCellColor(excelWorkbook, 'A1..BJ1', 6); % WT change cell color of first row to yellow
    excelWorkbook.Save;                                       % Save the current state of the workbook.
    excelWorkbook.Close;                                      % Close the workbook.  Excel will stay open but be hidden.
    excelObject.Quit;                                         % Shut down the Excel server instance.
    clear('excelObject', 'excelWorkbook', 'yourFormula');     % Clear the excel object variable from MATLAB's memory.
    fprintf('Done interacting with Excel.\n');                % The clear finally shuts down the server and it no longer appears in Task Manager.
    

end



end



% --------------------------------------------------------------------
% 04 country-level EMF - young2ref
% --------------------------------------------------------------------


if flags.valSTC_04 == 1


% vector
table_15  = table(round(squeeze(EMF_15( 1, :, ages==age_young))', 1),...
                  round(squeeze(EMF_15( 2, :, ages==age_young))', 1),...
                  round(squeeze(EMF_15( 3, :, ages==age_young))', 1),...
                  round(squeeze(EMF_15( 4, :, ages==age_young))', 1),...
                  round(squeeze(EMF_15( 5, :, ages==age_young))', 1),...
                  round(squeeze(EMF_15( 6, :, ages==age_young))', 1),...
                  round(squeeze(EMF_15( 7, :, ages==age_young))', 1),...
                  'RowNames',countries.name, 'VariableNames', extremes_legend);
table_20  = table(round(squeeze(EMF_20( 1, :, ages==age_young))', 1),...
                  round(squeeze(EMF_20( 2, :, ages==age_young))', 1),...
                  round(squeeze(EMF_20( 3, :, ages==age_young))', 1),...
                  round(squeeze(EMF_20( 4, :, ages==age_young))', 1),...
                  round(squeeze(EMF_20( 5, :, ages==age_young))', 1),...
                  round(squeeze(EMF_20( 6, :, ages==age_young))', 1),...
                  round(squeeze(EMF_20( 7, :, ages==age_young))', 1),...
                  'RowNames',countries.name, 'VariableNames', extremes_legend);
table_NDC = table(round(squeeze(EMF_NDC(1, :, ages==age_young))', 1),...
                  round(squeeze(EMF_NDC(2, :, ages==age_young))', 1),...
                  round(squeeze(EMF_NDC(3, :, ages==age_young))', 1),...
                  round(squeeze(EMF_NDC(4, :, ages==age_young))', 1),...
                  round(squeeze(EMF_NDC(5, :, ages==age_young))', 1),...
                  round(squeeze(EMF_NDC(6, :, ages==age_young))', 1),...
                  round(squeeze(EMF_NDC(7, :, ages==age_young))', 1),...
                  'RowNames',countries.name, 'VariableNames', extremes_legend);


% vector - save as txt file
writetable(table_15 , 'follow-up_STC\data\04_EMF_young2ref_15.xlsx' , 'WriteRowNames', true);
writetable(table_20 , 'follow-up_STC\data\04_EMF_young2ref_20.xlsx' , 'WriteRowNames', true);
writetable(table_NDC, 'follow-up_STC\data\04_EMF_young2ref_NDC.xlsx', 'WriteRowNames', true);


end



% --------------------------------------------------------------------
% 05 EMF young2ref excel sheets
% --------------------------------------------------------------------


if flags.valSTC_05 == 1   


    
% loop over GMT steps
for ind_BE=1:nGMTsteps

    
    % get multi-model mean EMF relative to same birth cohort living in *high-income country* under the same scenario
    EMF_perregion_lowincome2highincome_BE(:,:,:,ind_BE) = exposure_perregion_BE(:,:,:,ind_BE) ./ exposure_perregion_BE( :, 3, :, ind_BE);


end
    
    
% loop over extreme event categories
for ind_extreme = 1:nextremes+1
% for ind_extreme = nextremes+1


    % get associated Exposure Multiplication Factors (EMF) - LOW INCOME VS HIGH INCOME
    EMF_region_plot_BE = round( squeeze(EMF_perregion_lowincome2highincome_BE(ind_extreme, 5, :, :))', 1);


    % prepare output matrix
    matrix_out_05                     = NaN(nGMTsteps+1, nbirthyears+1);
    matrix_out_05(2:nGMTsteps+1, 1)   = fliplr(round(GMT_BE(end,:), 1));
    matrix_out_05(1, 2:nbirthyears+1) = ages;
    matrix_out_05(2:end, 2:end)       = flipud(EMF_region_plot_BE); 


    % save as excell file
    writematrix(matrix_out_05, 'follow-up_STC\data\05_EMF_lowincome2highincome.xlsx', 'Sheet', extremes_legend{ind_extreme});


end


% change cell color of first row and column to yellow
excelObject = actxserver('Excel.Application');            % Launch an Excel server using ActiveX (Windows ONLY).
fullFileName = fullfile(pwd, filename);                   % Create the filename of the existing workbook.
excelWorkbook = excelObject.workbooks.Open(fullFileName); % Open the workbook from disk.
excelObject.Visible = true;                               % Excel is invisible so far.  Make it visible.
Excel_utils.FormatCellColor(excelWorkbook, 'A1..A29', 6); % WT change cell color of first column to yellow
Excel_utils.FormatCellColor(excelWorkbook, 'A1..BJ1', 6); % WT change cell color of first row to yellow
excelWorkbook.Save;                                       % Save the current state of the workbook.
excelWorkbook.Close;                                      % Close the workbook.  Excel will stay open but be hidden.
excelObject.Quit;                                         % Shut down the Excel server instance.
clear('excelObject', 'excelWorkbook', 'yourFormula');     % Clear the excel object variable from MATLAB's memory.
fprintf('Done interacting with Excel.\n');                % The clear finally shuts down the server and it no longer appears in Task Manager.


end


% --------------------------------------------------------------------
% 06 global, all-extreme young2pic EMF (for COP26 exhibition)
% --------------------------------------------------------------------


if flags.valSTC_06 == 1


for ind_extreme = 1:nextremes+1
    
    % get matrix - 1.5
    valc_EMF_young2pic_15( :,1) = birth_years;
    valc_EMF_young2pic_15( :,2) = squeeze(round(EMF_perregion_15_young2pic( ind_extreme, 12, :), 2));


    % get matrix - NDC
    valc_EMF_young2pic_NDC(:,1) = birth_years;
    valc_EMF_young2pic_NDC(:,2) = squeeze(round(EMF_perregion_NDC_young2pic(ind_extreme, 12, :), 2));


    % vector - save as txt file
    writematrix(valc_EMF_young2pic_15 , ['follow-up_STC\data\06_EMF_young2pic_global_15_'  extremes_legend{ind_extreme} '.xlsx']);
    writematrix(valc_EMF_young2pic_NDC, ['follow-up_STC\data\06_EMF_young2pic_global_NDC_' extremes_legend{ind_extreme} '.xlsx']);

end

end



% --------------------------------------------------------------------
% 07 country-level EMF - young2pic, 
% --------------------------------------------------------------------


if flags.valSTC_07 == 1


% vector
table_15_young2pic  = table(round(squeeze(EMF_15_young2pic( 1, :, ages==age_young))', 1),...
                            round(squeeze(EMF_15_young2pic( 2, :, ages==age_young))', 1),...
                            round(squeeze(EMF_15_young2pic( 3, :, ages==age_young))', 1),...
                            round(squeeze(EMF_15_young2pic( 4, :, ages==age_young))', 1),...
                            round(squeeze(EMF_15_young2pic( 5, :, ages==age_young))', 1),...
                            round(squeeze(EMF_15_young2pic( 6, :, ages==age_young))', 1),...
                            round(squeeze(EMF_15_young2pic( 7, :, ages==age_young))', 1),...
                            'RowNames',countries.name, 'VariableNames', extremes_legend);
table_20_young2pic  = table(round(squeeze(EMF_20_young2pic( 1, :, ages==age_young))', 1),...
                            round(squeeze(EMF_20_young2pic( 2, :, ages==age_young))', 1),...
                            round(squeeze(EMF_20_young2pic( 3, :, ages==age_young))', 1),...
                            round(squeeze(EMF_20_young2pic( 4, :, ages==age_young))', 1),...
                            round(squeeze(EMF_20_young2pic( 5, :, ages==age_young))', 1),...
                            round(squeeze(EMF_20_young2pic( 6, :, ages==age_young))', 1),...
                            round(squeeze(EMF_20_young2pic( 7, :, ages==age_young))', 1),...
                            'RowNames',countries.name, 'VariableNames', extremes_legend);
table_NDC_young2pic = table(round(squeeze(EMF_NDC_young2pic(1, :, ages==age_young))', 1),...
                            round(squeeze(EMF_NDC_young2pic(2, :, ages==age_young))', 1),...
                            round(squeeze(EMF_NDC_young2pic(3, :, ages==age_young))', 1),...
                            round(squeeze(EMF_NDC_young2pic(4, :, ages==age_young))', 1),...
                            round(squeeze(EMF_NDC_young2pic(5, :, ages==age_young))', 1),...
                            round(squeeze(EMF_NDC_young2pic(6, :, ages==age_young))', 1),...
                            round(squeeze(EMF_NDC_young2pic(7, :, ages==age_young))', 1),...
                            'RowNames',countries.name, 'VariableNames', extremes_legend);


% vector - save as txt file
writetable(table_15_young2pic , 'follow-up_STC\data\07_EMF_young2pic_15.xlsx' , 'WriteRowNames', true);
writetable(table_20_young2pic , 'follow-up_STC\data\07_EMF_young2pic_20.xlsx' , 'WriteRowNames', true);
writetable(table_NDC_young2pic, 'follow-up_STC\data\07_EMF_young2pic_NDC.xlsx', 'WriteRowNames', true);


end



% --------------------------------------------------------------------
% 08 Portugal EMF - young2pic, different ages of plaintifs Duarte 
% (for Norwegian Human Rights institute NHRI)
% --------------------------------------------------------------------


% Define ages of plaintifs in the Duarte case in the year 2000
ages_Duarte = (21:-1:8)'; % for birth years 1999, 2000, 2003, 2004, 2008, 2012, their age on 31/12/2020


if flags.valNHRI_01 == 1


% vector
table_NDC_young2pic_Portugal = table(round(squeeze(EMF_NDC_young2pic(1, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...
                                     round(squeeze(EMF_NDC_young2pic(2, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...
                                     round(squeeze(EMF_NDC_young2pic(3, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...
                                     round(squeeze(EMF_NDC_young2pic(4, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...
                                     round(squeeze(EMF_NDC_young2pic(5, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...
                                     round(squeeze(EMF_NDC_young2pic(6, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...
                                     'RowNames',string(ages_Duarte) , 'VariableNames', extremes_legend(1:6));
%                                      round(squeeze(EMF_NDC_young2pic(7, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...


% vector - save as txt file
writetable(table_NDC_young2pic_Portugal, 'follow-up_NHRI\data\01_EMF_young2pic_NDC_Portugal.xlsx', 'WriteRowNames', true);


end



% --------------------------------------------------------------------
% 09 Portugal reduced exposure - young2pic, different ages of plaintifs
% Duarte (for Norwegian Human Rights institute NHRI)
% --------------------------------------------------------------------


if flags.valNHRI_02 == 1

    
% compute reduced exposure - note that we define this differently than reduction in intergenerational burden!!
valp_reduced_exposure_NDCto15 = (exposure_15 - exposure_NDC) ./ exposure_NDC .* 100;
valp_reduced_exposure_NDCto20 = (exposure_20 - exposure_NDC) ./ exposure_NDC .* 100;


% vector NDC to 1.5°C
table_reduced_exposure_NDCto15_Portugal = table(round(squeeze(valp_reduced_exposure_NDCto15(1, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto15(2, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto15(3, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto15(4, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto15(5, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto15(6, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          'RowNames',string(ages_Duarte) , 'VariableNames', extremes_legend(1:6));
%                                           round(squeeze(EMF_NDC_young2pic(7, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...


% vector NDC to 2.0°C
table_reduced_exposure_NDCto20_Portugal = table(round(squeeze(valp_reduced_exposure_NDCto20(1, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto20(2, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto20(3, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto20(4, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto20(5, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          round(squeeze(valp_reduced_exposure_NDCto20(6, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 0),...
                                          'RowNames',string(ages_Duarte) , 'VariableNames', extremes_legend(1:6));
%                                           round(squeeze(EMF_NDC_young2pic(7, strcmp(countries.name, 'Portugal'), ismember(ages, ages_Duarte))), 1),...


% vector - save as txt file
writetable(table_reduced_exposure_NDCto15_Portugal, 'follow-up_NHRI\data\02_reduced_exposure_NDCto15_Portugal.xlsx', 'WriteRowNames', true);
writetable(table_reduced_exposure_NDCto20_Portugal, 'follow-up_NHRI\data\02_reduced_exposure_NDCto20_Portugal.xlsx', 'WriteRowNames', true);
    


end



% --------------------------------------------------------------------
% 08 Switzerland EMF - young2pic, different ages of plaintifs Duarte 
% (for Veruska Muccione)
% --------------------------------------------------------------------


if flags.valVerus_01 == 1


% vector
table_NDC_young2pic_Switzerland = table(round(squeeze(EMF_NDC_young2pic(1, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        round(squeeze(EMF_NDC_young2pic(2, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        round(squeeze(EMF_NDC_young2pic(3, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        round(squeeze(EMF_NDC_young2pic(4, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        round(squeeze(EMF_NDC_young2pic(5, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        round(squeeze(EMF_NDC_young2pic(6, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        round(squeeze(EMF_NDC_young2pic(7, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                        'RowNames',string(ages) , 'VariableNames', extremes_legend);


% vector - save as txt file
writetable(table_NDC_young2pic_Switzerland, 'follow-up_Veruska\data\01_EMF_young2pic_NDC_Switzerland.xlsx', 'WriteRowNames', true);


end



% --------------------------------------------------------------------
% 09 Portugal reduced exposure - young2pic, different ages of plaintifs
% Duarte (for Norwegian Human Rights institute NHRI)
% --------------------------------------------------------------------


if flags.valVerus_02 == 1

    
% compute reduced exposure - note that we define this differently than reduction in intergenerational burden!!
valp_reduced_exposure_NDCto15 = (exposure_15 - exposure_NDC) ./ exposure_NDC .* 100;
valp_reduced_exposure_NDCto20 = (exposure_20 - exposure_NDC) ./ exposure_NDC .* 100;


% vector NDC to 1.5°C
table_reduced_exposure_NDCto15_Switzerland = table(round(squeeze(valp_reduced_exposure_NDCto15(1, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto15(2, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto15(3, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto15(4, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto15(5, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto15(6, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto15(7, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                                   'RowNames',string(ages) , 'VariableNames', extremes_legend);


% vector NDC to 2.0°C
table_reduced_exposure_NDCto20_Switzerland = table(round(squeeze(valp_reduced_exposure_NDCto20(1, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto20(2, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto20(3, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto20(4, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto20(5, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto20(6, strcmp(countries.name, 'Switzerland'), :)), 0),...
                                                   round(squeeze(valp_reduced_exposure_NDCto20(7, strcmp(countries.name, 'Switzerland'), :)), 1),...
                                                   'RowNames',string(ages) , 'VariableNames', extremes_legend);


% vector - save as txt file
writetable(table_reduced_exposure_NDCto15_Switzerland, 'follow-up_Veruska\data\02_reduced_exposure_NDCto15_Switzerland.xlsx', 'WriteRowNames', true);
writetable(table_reduced_exposure_NDCto20_Switzerland, 'follow-up_Veruska\data\02_reduced_exposure_NDCto20_Switzerland.xlsx', 'WriteRowNames', true);
    


end



% --------------------------------------------------------------------
% 01 data of figure 1 for carbon brief
% --------------------------------------------------------------------


if flags.valCB_01 == 1


if flags.plot_fig1 == 1

% panel a
figure_01_a_means(1,:) = years_SR15(1:end-4);
figure_01_a_means(2,:) = landfrac_15_plot(1:end-4);
figure_01_a_means(3,:) = landfrac_20_plot(1:end-4);
figure_01_a_means(4,:) = landfrac_NDC_plot(1:end-4);

figure_01_a_std(  1,:) = years_SR15(1:end-4);
figure_01_a_std(  2,:) = landfrac_15_mms_plot(1:end-4)./2;
figure_01_a_std(  3,:) = landfrac_20_mms_plot(1:end-4)./2;
figure_01_a_std(  4,:) = landfrac_NDC_mms_plot(1:end-4)./2;


% panel b
figure_01_b_abs      = exposure_bars;

figure_01_b_EMF(1,1) = round(EMF_plot_15( ages==age_ref_plot));
figure_01_b_EMF(1,2) = round(EMF_plot_20( ages==age_ref_plot));
figure_01_b_EMF(1,3) = round(EMF_plot_NDC(ages==age_ref_plot));
figure_01_b_EMF(2,1) = round(EMF_plot_15( ages==age_young));
figure_01_b_EMF(2,2) = round(EMF_plot_20( ages==age_young));
figure_01_b_EMF(2,3) = round(EMF_plot_NDC(ages==age_young));


% panel c
figure_01_c(1,:) = ages;
figure_01_c(2,:) = EMF_plot_15;
figure_01_c(3,:) = EMF_plot_20;
figure_01_c(4,:) = EMF_plot_NDC;


% save as txt file
writematrix(figure_01_a_means, 'follow-up_press\carbonbrief\figure_01_a_means.txt', 'Delimiter', 'tab');
writematrix(figure_01_a_std  , 'follow-up_press\carbonbrief\figure_01_a_std.txt'  , 'Delimiter', 'tab');
writematrix(figure_01_b_abs  , 'follow-up_press\carbonbrief\figure_01_b_abs.txt'  , 'Delimiter', 'tab');
writematrix(figure_01_b_EMF  , 'follow-up_press\carbonbrief\figure_01_b_EMF.txt'  , 'Delimiter', 'tab');
writematrix(figure_01_c      , 'follow-up_press\carbonbrief\figure_01_c.txt'      , 'Delimiter', 'tab');



end



end



% --------------------------------------------------------------------
% 02 data of supplementary figure 6 and 7 for carbon brief
% --------------------------------------------------------------------


if flags.valCB_02 == 1


if flags.plot_sfig6 == 1

    
% Raster - save lat and lon as txt file
writematrix(lat_mod, 'follow-up_press\carbonbrief\figure_S06_raster_lat.txt', 'Delimiter', 'tab');
writematrix(lon_mod, 'follow-up_press\carbonbrief\figure_S06_raster_lon.txt', 'Delimiter', 'tab');
    
    
% loop over extremes
for i=1:nextremes    
    
    % Raster - save as txt file - per extreme
    writematrix(EMF_map_young2ref_15(i).perextreme_plot , ['follow-up_press\carbonbrief\figure_S06_raster_15_'  extremes_legend{i} '.txt'], 'Delimiter', 'tab');
    writematrix(EMF_map_young2ref_NDC(i).perextreme_plot, ['follow-up_press\carbonbrief\figure_S06_raster_NDC_' extremes_legend{i} '.txt'], 'Delimiter', 'tab');
    
    
end


% Raster - save as txt file - all extremes
writematrix(EMF_map_young2ref_15_allhazards_plot , ['follow-up_press\carbonbrief\figure_S06_raster_15_'  extremes_legend{nextremes + 1} '.txt'], 'Delimiter', 'tab');
writematrix(EMF_map_young2ref_NDC_allhazards_plot, ['follow-up_press\carbonbrief\figure_S06_raster_NDC_' extremes_legend{nextremes + 1} '.txt'], 'Delimiter', 'tab');


% vector
table_15  = table(squeeze(EMF_15( 1, :, ages==age_young))',...
                  squeeze(EMF_15( 2, :, ages==age_young))',...
                  squeeze(EMF_15( 3, :, ages==age_young))',...
                  squeeze(EMF_15( 4, :, ages==age_young))',...
                  squeeze(EMF_15( 5, :, ages==age_young))',...
                  squeeze(EMF_15( 6, :, ages==age_young))',...
                  squeeze(EMF_15( 7, :, ages==age_young))',...
                  'RowNames',countries.name, 'VariableNames', extremes_legend);
table_NDC = table(squeeze(EMF_NDC(1, :, ages==age_young))',...
                  squeeze(EMF_NDC(2, :, ages==age_young))',...
                  squeeze(EMF_NDC(3, :, ages==age_young))',...
                  squeeze(EMF_NDC(4, :, ages==age_young))',...
                  squeeze(EMF_NDC(5, :, ages==age_young))',...
                  squeeze(EMF_NDC(6, :, ages==age_young))',...
                  squeeze(EMF_NDC(7, :, ages==age_young))',...
                  'RowNames',countries.name, 'VariableNames', extremes_legend);


% vector - save as txt file
writetable(table_15 , 'follow-up_press\carbonbrief\figure_S06_vector_15.txt' , 'Delimiter', 'tab', 'WriteRowNames', true);
writetable(table_NDC, 'follow-up_press\carbonbrief\figure_S06_vector_NDC.txt', 'Delimiter', 'tab', 'WriteRowNames', true);


end



end



% --------------------------------------------------------------------
% 03 data of supplementary figure 9 and 10 for carbon brief
% --------------------------------------------------------------------


if flags.valCB_03 == 1


if flags.plot_sfig9 == 1

    
% store cohort size
cohort_size_rel = cell2mat(regions.cohort_size_rel);
table_cohort_size  = table(cohort_size_rel(:,1),...
                           cohort_size_rel(:,2),...
                           cohort_size_rel(:,3),...
                           cohort_size_rel(:,4),...
                           'RowNames', regions.name, 'VariableNames', {'60', '40', '20', '0'});

% save as txt file
writetable(table_cohort_size , 'follow-up_press\carbonbrief\figure_S09_cohort_size.txt' , 'Delimiter', 'tab', 'WriteRowNames', true);

                       
% loop over extremes
for i=1:nextremes    
    
    % Raster - save as txt file - per extreme
    writematrix(squeeze(EMF_perregion_15(i,:,:)) , ['follow-up_press\carbonbrief\figure_S09_15_'  extremes_legend{i} '.txt'], 'Delimiter', 'tab');
    writematrix(squeeze(EMF_perregion_20(i,:,:)) , ['follow-up_press\carbonbrief\figure_S09_20_'  extremes_legend{i} '.txt'], 'Delimiter', 'tab');
    writematrix(squeeze(EMF_perregion_NDC(i,:,:)), ['follow-up_press\carbonbrief\figure_S09_NDC_' extremes_legend{i} '.txt'], 'Delimiter', 'tab');
    
    
end


end


end



% --------------------------------------------------------------------
% IPCC numbers for Europe for IPCC AR6 (Europe Chapter)
% --------------------------------------------------------------------


if flags.valAR6 == 1


% get extra exposure for newborns under 3.5°C compared to 1.5°C
[~,indmin_15] =   min(abs(GMT_BE(end,:) - 1.5))
[~,indmin_35] =   min(abs(GMT_BE(end,:) - 3.5))
valc_EMF_EUCA_35to15_BE_allextremes = round(exposure_perregion_BE(7,2,ages==age_young, indmin_35) ./ exposure_perregion_BE(7,2, ages==age_young, indmin_15),1);
valc_EMF_EUCA_35to15_BE_perextreme = round(exposure_perregion_BE(:,2,ages==age_young, indmin_35) ./ exposure_perregion_BE(:,2, ages==age_young, indmin_15),2);


valc_EMF_EUCA_young2pic_35_allextremes = round(EMF_perregion_young2pic_BE(:,2,ages==age_young,indmin_35),2);


end



% --------------------------------------------------------------------
% BBC (Jessica Furst)
% --------------------------------------------------------------------


if flags.valBBC == 1


% get extra exposure for newborns under NDC in low and high-income countries
valc_EMF_NDC_allextremes_lowincome  = round(EMF_perregion_NDC(7,5,ages==age_young),1)
valc_EMF_NDC_allextremes_highincome = round(EMF_perregion_NDC(7,3,ages==age_young),1)


% Change in exposure per region (tropical cyclones) 
valc_EMF_young2ref_perregion_NDC_allextremes      = regions.name;
valc_EMF_young2ref_perregion_NDC_allextremes(:,2) = num2cell(round(squeeze(EMF_perregion_NDC(nextremes+1, :, ages == 0)), 1))


end



% --------------------------------------------------------------------
% The Guardian (Damian Carrington)
% --------------------------------------------------------------------


if flags.valTG == 1 && flags.plot_fig1 == 1


% numbers for fig 1b
valc_fig1a_lines  = [years_SR15(1:end-4) landfrac_15_plot(1:end-4)        landfrac_20_plot(1:end-4)         landfrac_NDC_plot(1:end-4)       ]
valc_fig1a_std    = [years_SR15(1:end-4) landfrac_15_mms_plot(1:end-4)./2 landfrac_20_mms_plot(1:end-4)./2  landfrac_NDC_mms_plot(1:end-4)./2]


% numbers for fig 1b
valc_fig1b  = exposure_bars


% numbers for fig 1c
valc_fig1c_lines  = [ages EMF_plot_NDC EMF_plot_20 EMF_plot_15]
valc_fig1c_bars   = [[EMF_plot_q25_NDC(end); EMF_plot_NDC(end); EMF_plot_q75_NDC(end)] ...
                     [EMF_plot_q25_20(end);  EMF_plot_20(end);  EMF_plot_q75_20(end)]  ...
                     [EMF_plot_q25_15(end);  EMF_plot_15(end);  EMF_plot_q75_15(end)  ]]

                 
% save data
writematrix(valc_fig1a_lines, 'follow-up_press\guardian\fig1a_lines.xlsx');
writematrix(valc_fig1a_std  , 'follow-up_press\guardian\fig1a_std.xlsx'  );
writematrix(valc_fig1b      , 'follow-up_press\guardian\fig1b.xlsx'      );
writematrix(valc_fig1c_lines, 'follow-up_press\guardian\fig1c_lines.xlsx');
writematrix(valc_fig1c_bars , 'follow-up_press\guardian\fig1c_bars.xlsx' );


end



% --------------------------------------------------------------------
% DPA (Gregor Bauernfeind)
% --------------------------------------------------------------------


if flags.valDPA == 1 


% EMF Germany young2ref - NDC
valc_EMF_young2ref_perextreme_NDC_Germany      = extremes_legend';
valc_EMF_young2ref_perextreme_NDC_Germany(:,2) = num2cell(round( EMF_NDC(:, find(contains(countries.name,'Germany')), ages==age_young)   , 2 )) % EMF


% EMF Germany young2ref - 1.5
valc_EMF_young2ref_perextreme_15_Germany      = extremes_legend';
valc_EMF_young2ref_perextreme_15_Germany(:,2) = num2cell(round( EMF_15( :, find(contains(countries.name,'Germany')), ages==age_young)   , 2 )) % EMF


% save data
writecell(valc_EMF_young2ref_perextreme_NDC_Germany, 'follow-up_press\dpa\EMF_young2ref_perextreme_NDC_Germany.xlsx');
writecell(valc_EMF_young2ref_perextreme_15_Germany , 'follow-up_press\dpa\EMF_young2ref_perextreme_15_Germany.xlsx' );


end



% --------------------------------------------------------------------
% Washington Post (Sarah Kaplan)
% --------------------------------------------------------------------


if flags.valWP == 1 


% 6-yr old at 3°C, youngref: 
valc_EMF_6yr_3deg_young2ref      = extremes_legend';
valc_EMF_6yr_3deg_young2ref(:,2) =  num2cell(round(squeeze(EMF_perregion_young2ref_BE(:, 12, ages == 6, GMT_steps == 3)), 1))



% aggregated exposure change of age cohorts below 20 under 1.5°
valc_EMF_allextremes_young2ref      = ages(ages<=20);
valc_EMF_allextremes_young2ref(:,2) = squeeze(EMF_perregion_young2ref_BE(nextremes + 1, 12, ages<=20, GMT_steps == 1.5));  % 1.5°C of warming
valc_EMF_allextremes_young2ref(:,3) = squeeze(EMF_perregion_young2ref_BE(nextremes + 1, 12, ages<=20, GMT_steps == 3))     % 3°C of warming


end



% --------------------------------------------------------------------
% ISIpedia (Julia Reimann, Martin Park, Mahé Perrette)
% --------------------------------------------------------------------


if flags.valISIpedia == 1  


% 01: raster and vectorised data
if flags.plot_sfig7 == 1
    
save('follow-up_press\isipedia\01_EMF_map_young2ref_raster', 'EMF_map_young2ref_15_allhazards'     , 'EMF_map_young2ref_20_allhazards'     , 'EMF_map_young2ref_NDC_allhazards'      ...
                                                           , 'EMF_map_young2ref_15_allhazards_plot', 'EMF_map_young2ref_20_allhazards_plot', 'EMF_map_young2ref_NDC_allhazards_plot' ...
                                                           , 'EMF_map_young2ref_15'                , 'EMF_map_young2ref_20'                , 'EMF_map_young2ref_NDC'                 ...
                                                           , 'lat_mod', 'lon_mod');
save('follow-up_press\isipedia\01_EMF_map_young2pic_raster', 'EMF_map_young2pic_15_allhazards'     , 'EMF_map_young2pic_20_allhazards'     , 'EMF_map_young2pic_NDC_allhazards'      ...
                                                           , 'EMF_map_young2pic_15_allhazards_plot', 'EMF_map_young2pic_20_allhazards_plot', 'EMF_map_young2pic_NDC_allhazards_plot' ...
                                                           , 'EMF_map_young2pic_15'                , 'EMF_map_young2pic_20'                , 'EMF_map_young2pic_NDC'                 ...
                                                           , 'lat_mod', 'lon_mod');
save('follow-up_press\isipedia\01_exposure_map_raster'     , 'exposure_map_15_allhazards'          , 'exposure_map_20_allhazards'          , 'exposure_map_NDC_allhazards'           ...
                                                           , 'exposure_map_15_allhazards_plot'     , 'exposure_map_20_allhazards_plot'     , 'exposure_map_NDC_allhazards_plot'      ...
                                                           , 'exposure_map_15'                     , 'exposure_map_20'                     , 'exposure_map_NDC'                      ...
                                                           , 'lat_mod', 'lon_mod');
                                                       
save('follow-up_press\isipedia\01_EMF_map_young2ref_vector', 'ages', 'birth_years', 'extremes_legend', 'countries', 'EMF_15'          , 'EMF_20'          , 'EMF_NDC'          );
save('follow-up_press\isipedia\01_EMF_map_young2pic_vector', 'ages', 'birth_years', 'extremes_legend', 'countries', 'EMF_15_young2pic', 'EMF_20_young2pic', 'EMF_NDC_young2pic');
save('follow-up_press\isipedia\01_exposure_map_vector'     , 'ages', 'birth_years', 'extremes_legend', 'countries', 'exposure_15'     , 'exposure_20'     , 'exposure_NDC'     );

end





% 02: burning embers data
save('follow-up_press\isipedia\02_EMF_burning_embers', 'ages', 'birth_years', 'extremes_legend', 'regions', 'GMT_BE', 'EMF_perregion_young2pic_BE');


% 03: line plots per region
save('follow-up_press\isipedia\03_EMF_line_perregion_young2pic', 'ages', 'birth_years', 'extremes_legend', 'regions', 'EMF_perregion_15_young2pic', 'EMF_perregion_20_young2pic', 'EMF_perregion_NDC_young2pic');
save('follow-up_press\isipedia\03_EMF_line_perregion_young2ref', 'ages', 'birth_years', 'extremes_legend', 'regions', 'EMF_perregion_15'          , 'EMF_perregion_20'          , 'EMF_perregion_NDC'          );


% 04: line plots per country
save('follow-up_press\isipedia\04_EMF_line_percountry_young2ref', 'ages', 'birth_years', 'extremes_legend', 'countries', 'EMF_15'          , 'EMF_20'          , 'EMF_NDC'          );




end



% --------------------------------------------------------------------
% My Climate Future (Jonas Parnow)
% --------------------------------------------------------------------


if flags.valMCF == 1  


% loop over regions
for ind_region = 1:12
% for ind_region = 12
  
    
    % loop over extreme event categories
    for ind_extreme = 1:nextremes
%     for ind_extreme = nextremes

        
        % get associated Exposure Multiplication Factors (EMF)
        EMF_region_plot_BE = round(squeeze(EMF_perregion_young2pic_BE(ind_extreme, ind_region, :, :))', 1);

        
        % prepare output matrix
        matrix_out_02                     = NaN(nGMTsteps+1, nbirthyears+1);
        matrix_out_02(2:nGMTsteps+1, 1)   = fliplr(round(GMT_BE(end,:), 2));
        matrix_out_02(1, 2:nbirthyears+1) = ages;
        matrix_out_02(2:end, 2:end)       = flipud(EMF_region_plot_BE); 
        
        
        % save as excell file
        writematrix(matrix_out_02, ['follow-up_myclimatefuture\data\02_EMF_young2pic_' extremes{ind_extreme} '_' regions.name{ind_region} '.csv']);

        
    end


end


% get list of countries in each region (for the list on the website)
for ind_region = 1:12

    % get number of countries per region
    nmembercountries = length(regions.member_countries{ind_region});
   
    % loop over member countries and store them in one string
    membercountrystring{ind_region,1} = cell2mat(regions.member_countries{ind_region}(1));
    for ind_membercountry = 2:nmembercountries
        membercountrystring{ind_region,1} = [membercountrystring{ind_region,1} ', ' cell2mat(regions.member_countries{ind_region}(ind_membercountry))];
    end
        
end
membercountrystring = membercountrystring


end


% --------------------------------------------------------------------
% Scientific american EMF young2pic
% --------------------------------------------------------------------


if flags.valsciam == 1


% loop over regions
% for ind_region = 1:12
for ind_region = 12
  
    
    % loop over extreme event categories
    for ind_extreme = 1:nextremes+1
%     for ind_extreme = nextremes+1

        
        % get associated Exposure Multiplication Factors (EMF)
        EMF_region_plot_BE = round(squeeze(EMF_perregion_young2pic_BE(ind_extreme, ind_region, :, :))', 1);

        
        % prepare output matrix
        valc_EMF_young2pic_global_perextreme                     = NaN(nGMTsteps+1, nbirthyears+1);
        valc_EMF_young2pic_global_perextreme(2:nGMTsteps+1, 1)   = fliplr(round(GMT_BE(end,:), 1));
        valc_EMF_young2pic_global_perextreme(1, 2:nbirthyears+1) = ages;
        valc_EMF_young2pic_global_perextreme(2:end, 2:end)       = flipud(EMF_region_plot_BE); 
        
        
        % save as excell file
        filename = ['follow-up_press\scientific_american\EMF_young2pic_' regions.name{ind_region} '.xlsx'];
        writematrix(valc_EMF_young2pic_global_perextreme, filename, 'Sheet', extremes_legend{ind_extreme});
 
        
    end

         
    % change cell color of first row and column to yellow
    excelObject = actxserver('Excel.Application');            % Launch an Excel server using ActiveX (Windows ONLY).
    fullFileName = fullfile(pwd, filename);                   % Create the filename of the existing workbook.
    excelWorkbook = excelObject.workbooks.Open(fullFileName); % Open the workbook from disk.
    excelObject.Visible = true;                               % Excel is invisible so far.  Make it visible.
    Excel_utils.FormatCellColor(excelWorkbook, 'A1..A29', 6); % WT change cell color of first column to yellow
    Excel_utils.FormatCellColor(excelWorkbook, 'A1..BJ1', 6); % WT change cell color of first row to yellow
    excelWorkbook.Save;                                       % Save the current state of the workbook.
    excelWorkbook.Close;                                      % Close the workbook.  Excel will stay open but be hidden.
    excelObject.Quit;                                         % Shut down the Excel server instance.
    clear('excelObject', 'excelWorkbook', 'yourFormula');     % Clear the excel object variable from MATLAB's memory.
    fprintf('Done interacting with Excel.\n');                % The clear finally shuts down the server and it no longer appears in Task Manager.
    

end



if flags.plot_sfig10 == 1
    
    % get relative cohort sizes
%     valc_cohort_size_rel = cell2mat(regions.cohort_size_rel);
    
    clear valc_cohort_size_rel
    valc_cohort_size_rel(1,2:5)            = num2cell(ages_pie);
    valc_cohort_size_rel(2:nregions+1,1)   = regions.name';
    valc_cohort_size_rel(2:nregions+1,2:5) = regions.cohort_size_rel;

    % save as excell file
    writecell(valc_cohort_size_rel, 'follow-up_press\scientific_american\cohort_size.xlsx');

    
    % loop over 0, 20, 40, 60
    for ind_age_pie=1:length(ages_pie)
        
        % store data in matrix and save it - NDC
        valc_EMF_perregion_NDC_young2pic_pies(1,2:nregions+1)             = regions.name;
        valc_EMF_perregion_NDC_young2pic_pies(2:nextremes+2,1)            = extremes_legend;
        valc_EMF_perregion_NDC_young2pic_pies(2:nextremes+2,2:nregions+1) = num2cell(EMF_perregion_NDC_young2pic(:,:,ages==ages_pie(ind_age_pie)));       
        writecell(valc_EMF_perregion_NDC_young2pic_pies, ['follow-up_press\scientific_american\EMF_perregion_NDC_young2pic_' num2str(ages_pie(ind_age_pie)) '.xlsx' ]);
        
        % store data in matrix and save it - 2.0°C
        valc_EMF_perregion_20_young2pic_pies(1,2:nregions+1)             = regions.name;
        valc_EMF_perregion_20_young2pic_pies(2:nextremes+2,1)            = extremes_legend;
        valc_EMF_perregion_20_young2pic_pies(2:nextremes+2,2:nregions+1) = num2cell(EMF_perregion_20_young2pic(:,:,ages==ages_pie(ind_age_pie)));       
        writecell(valc_EMF_perregion_20_young2pic_pies, ['follow-up_press\scientific_american\EMF_perregion_20_young2pic_' num2str(ages_pie(ind_age_pie)) '.xlsx' ]);
        
        % store data in matrix and save it - 1.5°C
        valc_EMF_perregion_15_young2pic_pies(1,2:nregions+1)             = regions.name;
        valc_EMF_perregion_15_young2pic_pies(2:nextremes+2,1)            = extremes_legend;
        valc_EMF_perregion_15_young2pic_pies(2:nextremes+2,2:nregions+1) = num2cell(EMF_perregion_15_young2pic(:,:,ages==ages_pie(ind_age_pie)));       
        writecell(valc_EMF_perregion_15_young2pic_pies, ['follow-up_press\scientific_american\EMF_perregion_15_young2pic_' num2str(ages_pie(ind_age_pie)) '.xlsx' ]);
    end
    
    
    
end


end




% --------------------------------------------------------------------
% Values for Bloom
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valBloom == 1


    % change in GMT from emissions 2.4-2.8 GtCO2e
    CO2_emissions_Bloom_min  = 2.4E9;   % taking the middle of the range
    CO2_emissions_Bloom_mean = 2.6E9;   % taking the middle of the range
    CO2_emissions_Bloom_max  = 2.8E9;   % taking the middle of the range
    dGMT_Bloom_min           = 0.5 ./ 1000E9 * CO2_emissions_Bloom_min    % '1000 GT = 0.5 °C rule'  ==> IMPROVE THIS PART
    dGMT_Bloom_mean          = 0.5 ./ 1000E9 * CO2_emissions_Bloom_mean   % '1000 GT = 0.5 °C rule'  ==> IMPROVE THIS PART
    dGMT_Bloom_max           = 0.5 ./ 1000E9 * CO2_emissions_Bloom_max    % '1000 GT = 0.5 °C rule'  ==> IMPROVE THIS PART
    
    
    % extract the lifetime absolute number of heatwaves that newborns will
    % experience under the burning embers pathways:
    valc_exposure_heatwaves_newborns     = squeeze(exposure_perregion_BE(5,12,ages==0,:));  % lifetime absolute heatwave exposure for newborns
    valc_GMT_2100                        = squeeze(GMT_BE(end,:));                          % 2100 GMT anomaly

    
    % fit a linear curve through the relation and get the slope
    valc_pf                       = polyfit(valc_GMT_2100,valc_exposure_heatwaves_newborns,1);
    valc_slope_exposure_heatwaves = valc_pf(1);


    % exctract the number of newborns
    nr_newborns                   = valp_cohort_size_abs(61,12);

    % get the diagnostics
    nr_extra_heatwaves_newborns_Bloom_min  = valc_slope_exposure_heatwaves * dGMT_Bloom_min ;
    nr_extra_heatwaves_newborns_Bloom_mean = valc_slope_exposure_heatwaves * dGMT_Bloom_mean;
    nr_extra_heatwaves_newborns_Bloom_max  = valc_slope_exposure_heatwaves * dGMT_Bloom_max ;
    nr_children_facing_extra_heatwave_Bloom_min   = floor(nr_newborns * nr_extra_heatwaves_newborns_Bloom_min )
    nr_children_facing_extra_heatwave_Bloom_mean  = floor(nr_newborns * nr_extra_heatwaves_newborns_Bloom_mean)
    nr_children_facing_extra_heatwave_Bloom_max   = floor(nr_newborns * nr_extra_heatwaves_newborns_Bloom_max )
   
    
    
    % second calculation: nr of heat-related deaths between today and 2100
    mortality_cost_carbon = 4434;                                         % “4,434 metric tons of carbon dioxide in 2020 […] causes one excess death globally in expectation between 2020-2100” (Bressler, 2021)
    mortality_Bloom_min   = CO2_emissions_Bloom_min  ./ mortality_cost_carbon  % we find 541270
    mortality_Bloom_mean  = CO2_emissions_Bloom_mean ./ mortality_cost_carbon  % we find 586380
    mortality_Bloom_max   = CO2_emissions_Bloom_max  ./ mortality_cost_carbon  % we find 631480

end



% --------------------------------------------------------------------
% Values for testimony in Norwegian court case
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valNorLic == 1


    % change in GMT from emissions of three oil/gas fields
    dGMT_Norway = 0.00023;
    
    
    % extract the lifetime absolute number of heatwaves that newborns will
    % experience under the burning embers pathways:
    valc_exposure_heatwaves_newborns     = squeeze(exposure_perregion_BE(5,12,ages==0,:));  % lifetime absolute heatwave exposure for newborns
    valc_GMT_2100                        = squeeze(GMT_BE(end,:));                          % 2100 GMT anomaly

    
    % fit a linear curve through the relation and get the slope
    valc_pf                       = polyfit(valc_GMT_2100,valc_exposure_heatwaves_newborns,1);
    valc_slope_exposure_heatwaves = valc_pf(1);


    % exctract the number of newborns
    nr_newborns                              = valp_cohort_size_abs(61,12)


    % get the diagnostics
    nr_extra_heatwaves_newborns_from_FFfield = valc_slope_exposure_heatwaves * dGMT_Norway;
    nr_children_facing_extra_heatwave_Norway = floor(nr_newborns * nr_extra_heatwaves_newborns_from_FFfield)
    
    
    % plot the results
    valc_polyval_exposure_heatwaves = polyval(valc_pf,valc_GMT_2100);
    close all
    figure
    plot(valc_GMT_2100,valc_exposure_heatwaves_newborns); hold on % modelled 
    plot(valc_GMT_2100,valc_polyval_exposure_heatwaves)                                         %

end





% --------------------------------------------------------------------
% Values for Norway -- injunction order
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valNorStop == 1

    
    % Compute change in GMT from annual emissions
    CO2_emissions_Breidablikk_2024 = 8.8E6;                                                           % 8,9 MtCO2e from Breidablikk in 2024
    CO2_emissions_Breidablikk_2025 = 9.1E6;                                                           % 9,1 MtCO2e from Breidablikk in 2025
    CO2_emissions_Breidablikk_tot  = CO2_emissions_Breidablikk_2024 + CO2_emissions_Breidablikk_2025; % sum of 2024 and 2025 (8,9 + 9,1)
    dGMT_Breidablikk_2024          = TCRE * CO2_emissions_Breidablikk_2024                            
    dGMT_Breidablikk_2025          = TCRE * CO2_emissions_Breidablikk_2025 
    dGMT_Breidablikk_tot           = TCRE * CO2_emissions_Breidablikk_tot  
    
    
    for i=1:6 % loop over birth years 2020 to 2015
    
        % extract the lifetime absolute number of heatwaves that newborns will
        % experience under the burning embers pathways:
        valc_exposure_heatwaves_newborns = squeeze(exposure_perregion_BE(5,12,ages==ages(end-i+1),:));  % lifetime absolute heatwave exposure for newborns
        valc_GMT_2100                    = squeeze(GMT_BE(end,:));                                      % 2100 GMT anomaly


        % fit a linear curve through the relation and get the slope
        valc_pf                       = polyfit(valc_GMT_2100,valc_exposure_heatwaves_newborns,1);
        valc_slope_exposure_heatwaves = valc_pf(1);


        % exctract the number of newborns
        nr_newborns(i,1) = valp_cohort_size_abs(end-i+1,12);


        % Compute the average change in lifetime heatwave exposure for each birth cohort member
        nr_extra_heatwaves_newborns_Breidablikk_2024(i,1) = valc_slope_exposure_heatwaves * dGMT_Breidablikk_2024;
        nr_extra_heatwaves_newborns_Breidablikk_2025(i,1) = valc_slope_exposure_heatwaves * dGMT_Breidablikk_2025;
        nr_extra_heatwaves_newborns_Breidablikk_tot(i,1)  = valc_slope_exposure_heatwaves * dGMT_Breidablikk_tot ;


        % from the previous number, compute the number of children that are expected to face just one extra heatwave
        nr_children_facing_extra_heatwave_Breidablikk_2024(i,1) = floor(nr_newborns(i,1) * nr_extra_heatwaves_newborns_Breidablikk_2024(i,1)); % we find 5474 for 2020 birth cohort
        nr_children_facing_extra_heatwave_Breidablikk_2025(i,1) = floor(nr_newborns(i,1) * nr_extra_heatwaves_newborns_Breidablikk_2025(i,1)); % we find 5597 for 2020 birth cohort
        nr_children_facing_extra_heatwave_Breidablikk_tot(i,1)  = floor(nr_newborns(i,1) * nr_extra_heatwaves_newborns_Breidablikk_tot(i,1) ); % we find 11071 for 2020 birth cohort
       
    end
    
    
    % get the sum across the years 2015-2020
    nr_children_facing_extra_heatwave_Breidablikk_2024(7,1) = sum(nr_children_facing_extra_heatwave_Breidablikk_2024(1:6,1))
    nr_children_facing_extra_heatwave_Breidablikk_2025(7,1) = sum(nr_children_facing_extra_heatwave_Breidablikk_2025(1:6,1))
    nr_children_facing_extra_heatwave_Breidablikk_tot( 7,1) = sum(nr_children_facing_extra_heatwave_Breidablikk_tot( 1:6,1))
    
    
    % second calculation: nr of heat-related deaths between today and 2100
    mortality_Breidablikk_2024 = floor(CO2_emissions_Breidablikk_2024 ./ mortality_cost_carbon) % we find 2007
    mortality_Breidablikk_2025 = floor(CO2_emissions_Breidablikk_2025 ./ mortality_cost_carbon) % we find 2052
    mortality_Breidablikk_tot  = floor(CO2_emissions_Breidablikk_tot  ./ mortality_cost_carbon) % we find 4059

end






% --------------------------------------------------------------------
% Values for Norway -- appeal phase
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valNorAppea == 1


    % TOTAL emissions (MtCO2e)
    CO2_emissions_Tyrving_tot     = 12E6; 
    CO2_emissions_Breidablikk_tot = 87E6; 
    CO2_emissions_Yggdrasil_tot   = 365E6; 
    CO2_emissions_combined_tot    = CO2_emissions_Tyrving_tot + CO2_emissions_Breidablikk_tot + CO2_emissions_Yggdrasil_tot; 

    
    % ANNUAL emissions (MtCO2e)
    CO2_emissions_Tyrving_2024     =  0.435E6;                                                           % 8,9 MtCO2e from Breidablikk in 2024
    CO2_emissions_Tyrving_2025     =  3.170E6;                                                           % 9,1 MtCO2e from Breidablikk in 2025
    CO2_emissions_Breidablikk_2024 = 10.143E6;                                                           % 8,9 MtCO2e from Breidablikk in 2024
    CO2_emissions_Breidablikk_2025 = 10.143E6;                                                           % 9,1 MtCO2e from Breidablikk in 2025
    CO2_emissions_Yggdrasil_2027   = 28.800E6;                                                           % 8,9 MtCO2e from Breidablikk in 2024
    CO2_emissions_Yggdrasil_2028   = 42.400E6;                                                           % 9,1 MtCO2e from Breidablikk in 2025

    
    % Question 1: call function to compute the number of people facing one
    % extra heatwave due to the total emissions of the three fields
    valc_nr_children_facing_extra_heatwave_tot(:,1)    = mf_emissions2npeople(CO2_emissions_Tyrving_tot    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_tot(:,2)    = mf_emissions2npeople(CO2_emissions_Breidablikk_tot, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_tot(:,3)    = mf_emissions2npeople(CO2_emissions_Yggdrasil_tot  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_tot(:,4)    = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5)  % 5 is for heatwaves
    

    % Question 2: call function to compute the number of people facing one
    % extra ... due to the total emissions of the three fields
    % Droughts
    valc_nr_children_facing_extra_drought_tot(:,1)     = mf_emissions2npeople(CO2_emissions_Tyrving_tot    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_tot(:,2)     = mf_emissions2npeople(CO2_emissions_Breidablikk_tot, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_tot(:,3)     = mf_emissions2npeople(CO2_emissions_Yggdrasil_tot  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_tot(:,4)     = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3)  % 3 is for Droughts

    % Crop Failures
    valc_nr_children_facing_extra_cropfailure_tot(:,1) = mf_emissions2npeople(CO2_emissions_Tyrving_tot    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_tot(:,2) = mf_emissions2npeople(CO2_emissions_Breidablikk_tot, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_tot(:,3) = mf_emissions2npeople(CO2_emissions_Yggdrasil_tot  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_tot(:,4) = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2)  % 2 is for Crop Failure

    % Wildfires
    valc_nr_children_facing_extra_wildfire_tot(:,1)    = mf_emissions2npeople(CO2_emissions_Tyrving_tot    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_tot(:,2)    = mf_emissions2npeople(CO2_emissions_Breidablikk_tot, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_tot(:,3)    = mf_emissions2npeople(CO2_emissions_Yggdrasil_tot  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_tot(:,4)    = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1)  % 1 is for wildfires

    % Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_tot(:,1) = mf_emissions2npeople(CO2_emissions_Tyrving_tot    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_tot(:,2) = mf_emissions2npeople(CO2_emissions_Breidablikk_tot, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_tot(:,3) = mf_emissions2npeople(CO2_emissions_Yggdrasil_tot  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_tot(:,4) = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6)  % 6 is for Tropical cyclones

    % River floods
    valc_nr_children_facing_extra_riverflood_tot(:,1)  = mf_emissions2npeople(CO2_emissions_Tyrving_tot    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_tot(:,2)  = mf_emissions2npeople(CO2_emissions_Breidablikk_tot, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_tot(:,3)  = mf_emissions2npeople(CO2_emissions_Yggdrasil_tot  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_tot(:,4)  = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4)  % 4 is for River floods
    
    
    
    
    % just for plotting purposes
    valc_nr_children_facing_extra_heatwave_ref(:,4)    = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 0, 5);  % 5 is for heatwaves
    valc_nr_children_facing_extra_drought_ref(:,4)     = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 0, 3);  % 3 is for Droughts
    valc_nr_children_facing_extra_cropfailure_ref(:,4) = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 0, 2);  % 2 is for Crop Failure
    valc_nr_children_facing_extra_wildfire_ref(:,4)    = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 0, 1);  % 1 is for wildfires
    valc_nr_children_facing_extra_tropcyclone_ref(:,4) = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 0, 6);  % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_riverflood_ref(:,4)  = mf_emissions2npeople(CO2_emissions_combined_tot   , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 0, 4);  % 4 is for River floods

    
    % Question 3: compute the number of heat-related deaths between today and 2100 due to the total emissions of the three fields
    valc_mortality_tot(1,1) = floor(CO2_emissions_Tyrving_tot     ./ mortality_cost_carbon); % we find 2706
    valc_mortality_tot(1,2) = floor(CO2_emissions_Breidablikk_tot ./ mortality_cost_carbon); % we find 19621
    valc_mortality_tot(1,3) = floor(CO2_emissions_Yggdrasil_tot   ./ mortality_cost_carbon); % we find 82318
    valc_mortality_tot(1,4) = floor(CO2_emissions_combined_tot    ./ mortality_cost_carbon)  % we find 104645

    
    % Question 4: repeat questions 1-4 but for the annual emissions of individual fields
    % heatwaves
    valc_nr_children_facing_extra_heatwave_annual(:,1)    = mf_emissions2npeople(CO2_emissions_Tyrving_2024    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_annual(:,2)    = mf_emissions2npeople(CO2_emissions_Tyrving_2025    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_annual(:,3)    = mf_emissions2npeople(CO2_emissions_Breidablikk_2024, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_annual(:,4)    = mf_emissions2npeople(CO2_emissions_Breidablikk_2025, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_annual(:,5)    = mf_emissions2npeople(CO2_emissions_Yggdrasil_2027  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_annual(:,6)    = mf_emissions2npeople(CO2_emissions_Yggdrasil_2028  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 5)  % 5 is for heatwaves
    
    % Droughts
    valc_nr_children_facing_extra_drought_annual(:,1)     = mf_emissions2npeople(CO2_emissions_Tyrving_2024    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_annual(:,2)     = mf_emissions2npeople(CO2_emissions_Tyrving_2025    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_annual(:,3)     = mf_emissions2npeople(CO2_emissions_Breidablikk_2024, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_annual(:,4)     = mf_emissions2npeople(CO2_emissions_Breidablikk_2025, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_annual(:,5)     = mf_emissions2npeople(CO2_emissions_Yggdrasil_2027  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_annual(:,6)     = mf_emissions2npeople(CO2_emissions_Yggdrasil_2028  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 3)  % 3 is for Droughts
    
    % Crop Failures
    valc_nr_children_facing_extra_cropfailure_annual(:,1) = mf_emissions2npeople(CO2_emissions_Tyrving_2024    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_annual(:,2) = mf_emissions2npeople(CO2_emissions_Tyrving_2025    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_annual(:,3) = mf_emissions2npeople(CO2_emissions_Breidablikk_2024, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_annual(:,4) = mf_emissions2npeople(CO2_emissions_Breidablikk_2025, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_annual(:,5) = mf_emissions2npeople(CO2_emissions_Yggdrasil_2027  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_annual(:,6) = mf_emissions2npeople(CO2_emissions_Yggdrasil_2028  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 2)  % 2 is for Crop Failure
    
    % Wildfires
    valc_nr_children_facing_extra_wildfire_annual(:,1)    = mf_emissions2npeople(CO2_emissions_Tyrving_2024    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_annual(:,2)    = mf_emissions2npeople(CO2_emissions_Tyrving_2025    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_annual(:,3)    = mf_emissions2npeople(CO2_emissions_Breidablikk_2024, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_annual(:,4)    = mf_emissions2npeople(CO2_emissions_Breidablikk_2025, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_annual(:,5)    = mf_emissions2npeople(CO2_emissions_Yggdrasil_2027  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_annual(:,6)    = mf_emissions2npeople(CO2_emissions_Yggdrasil_2028  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 1)  % 1 is for wildfires
    
    % Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_annual(:,1) = mf_emissions2npeople(CO2_emissions_Tyrving_2024    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_annual(:,2) = mf_emissions2npeople(CO2_emissions_Tyrving_2025    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_annual(:,3) = mf_emissions2npeople(CO2_emissions_Breidablikk_2024, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_annual(:,4) = mf_emissions2npeople(CO2_emissions_Breidablikk_2025, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_annual(:,5) = mf_emissions2npeople(CO2_emissions_Yggdrasil_2027  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_annual(:,6) = mf_emissions2npeople(CO2_emissions_Yggdrasil_2028  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 6)  % 6 is for Tropical cyclones
    
    % River floods
    valc_nr_children_facing_extra_riverflood_annual(:,1)  = mf_emissions2npeople(CO2_emissions_Tyrving_2024    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_annual(:,2)  = mf_emissions2npeople(CO2_emissions_Tyrving_2025    , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_annual(:,3)  = mf_emissions2npeople(CO2_emissions_Breidablikk_2024, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_annual(:,4)  = mf_emissions2npeople(CO2_emissions_Breidablikk_2025, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_annual(:,5)  = mf_emissions2npeople(CO2_emissions_Yggdrasil_2027  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_annual(:,6)  = mf_emissions2npeople(CO2_emissions_Yggdrasil_2028  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 0, 4)  % 4 is for River floods
    
    
    % compute the number of heat-related deaths between today and 2100 due to the annual emissions of the three fields
    valc_mortality_annual(1,1) = floor(CO2_emissions_Tyrving_2024     ./ mortality_cost_carbon); % we find 98
    valc_mortality_annual(1,2) = floor(CO2_emissions_Tyrving_2025     ./ mortality_cost_carbon); % we find 714
    valc_mortality_annual(1,3) = floor(CO2_emissions_Breidablikk_2024 ./ mortality_cost_carbon); % we find 2287
    valc_mortality_annual(1,4) = floor(CO2_emissions_Breidablikk_2025 ./ mortality_cost_carbon); % we find 2287
    valc_mortality_annual(1,5) = floor(CO2_emissions_Yggdrasil_2027   ./ mortality_cost_carbon); % we find 6495
    valc_mortality_annual(1,6) = floor(CO2_emissions_Yggdrasil_2028   ./ mortality_cost_carbon)  % we find 9562

    
    
end




% --------------------------------------------------------------------
% Values for Norway -- ECtHR
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valNorECtHR == 1


    % TOTAL emissions Barents Sea South (MtCO2e: add E6 to express as tCO2e)
    CO2_emissions_BarentsSea_South_min          = 2880E6; 
    CO2_emissions_BarentsSea_South_median       = 5184E6; 
    CO2_emissions_BarentsSea_South_max          = 6336E6; 

 
    % TOTAL emissions Barents Sea South East (MtCO2e: add E6 to express as tCO2e)
    CO2_emissions_BarentsSea_Southeast_min      = 132E6; 
    CO2_emissions_BarentsSea_Southeast_median   = 722E6; 
    CO2_emissions_BarentsSea_Southeast_max      = 1627E6; 

 
    % TOTAL emissions Barents Sea South (MtCO2e: add E6 to express as tCO2e)
    CO2_emissions_BarentsSea_Southeast_lowscen  = 106.9E6; 
    CO2_emissions_BarentsSea_Southeast_highscen = 388.0E6; 

 
    
    % Question 1: call function to compute the number of people facing one
    % extra heatwave due to the total emissions of the three fields
    clear valc_nr_children_facing_extra_heatwave_BarentsSea
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,1)    = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,2)    = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,3)    = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,4)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,5)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,6)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,7)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea(:,8)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    
    
    % include a small add-on table for the reference years 1960-1970
    clear valc_nr_children_facing_extra_heatwave_BarentsSea_ref
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref      = valc_nr_children_facing_extra_heatwave_BarentsSea_ref(end,:);
    valc_nr_children_facing_extra_heatwave_BarentsSea_ref(2,:) = round(valc_nr_children_facing_extra_heatwave_BarentsSea(end,:) ./ valc_nr_children_facing_extra_heatwave_BarentsSea_ref(1,:) .* 100);
        

    % Question 2: call function to compute the number of people facing one
    % extra ... due to the total emissions of the three fields
    % Droughts
    clear valc_nr_children_facing_extra_drought_BarentsSea
    valc_nr_children_facing_extra_drought_BarentsSea(:,1)     = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,2)     = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,3)     = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,4)     = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,5)     = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,6)     = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,7)     = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea(:,8)     = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts

    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_BarentsSea
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure

    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_BarentsSea
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,1)    = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,2)    = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,3)    = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,4)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,5)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,6)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,7)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea(:,8)    = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires

    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_BarentsSea
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones

    % River floods
    clear valc_nr_children_facing_extra_riverflood_BarentsSea
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,1)  = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,2)  = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,3)  = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,4)  = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,5)  = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,6)  = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,7)  = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea(:,8)  = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods

        
    % include a small add-on table for the reference years 1960-1970
    % Droughts
    clear valc_nr_children_facing_extra_drought_BarentsSea_ref
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_BarentsSea_ref      = valc_nr_children_facing_extra_drought_BarentsSea_ref(end,:);
    valc_nr_children_facing_extra_drought_BarentsSea_ref(2,:) = round(valc_nr_children_facing_extra_drought_BarentsSea(end,:) ./ valc_nr_children_facing_extra_drought_BarentsSea_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_BarentsSea_ref
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref      = valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(end,:);
    valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(2,:) = round(valc_nr_children_facing_extra_cropfailure_BarentsSea(end,:) ./ valc_nr_children_facing_extra_cropfailure_BarentsSea_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_BarentsSea_ref
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref      = valc_nr_children_facing_extra_wildfire_BarentsSea_ref(end,:);
    valc_nr_children_facing_extra_wildfire_BarentsSea_ref(2,:) = round(valc_nr_children_facing_extra_wildfire_BarentsSea(end,:) ./ valc_nr_children_facing_extra_wildfire_BarentsSea_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref      = valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(end,:);
    valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(2,:) = round(valc_nr_children_facing_extra_tropcyclone_BarentsSea(end,:) ./ valc_nr_children_facing_extra_tropcyclone_BarentsSea_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % River floods
    clear valc_nr_children_facing_extra_riverflood_BarentsSea_ref
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,1) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_min         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,2) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_median      , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,3) = mf_emissions2npeople(CO2_emissions_BarentsSea_South_max         , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,4) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_min     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,5) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_median  , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,6) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_max     , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,7) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_lowscen , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(:,8) = mf_emissions2npeople(CO2_emissions_BarentsSea_Southeast_highscen, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref      = valc_nr_children_facing_extra_riverflood_BarentsSea_ref(end,:);
    valc_nr_children_facing_extra_riverflood_BarentsSea_ref(2,:) = round(valc_nr_children_facing_extra_riverflood_BarentsSea(end,:) ./ valc_nr_children_facing_extra_riverflood_BarentsSea_ref(1,:) .* 100);
        
    
    
    % Question 3: compute the number of heat-related deaths between today and 2100 due to the total emissions of the three fields
    valc_mortality_BarentsSea(1,1) = floor(CO2_emissions_BarentsSea_South_min          ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,2) = floor(CO2_emissions_BarentsSea_South_median       ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,3) = floor(CO2_emissions_BarentsSea_South_max          ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,4) = floor(CO2_emissions_BarentsSea_Southeast_min      ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,5) = floor(CO2_emissions_BarentsSea_Southeast_median   ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,6) = floor(CO2_emissions_BarentsSea_Southeast_max      ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,7) = floor(CO2_emissions_BarentsSea_Southeast_lowscen  ./ mortality_cost_carbon ./ 1000) .* 1000; % we find 
    valc_mortality_BarentsSea(1,8) = floor(CO2_emissions_BarentsSea_Southeast_highscen ./ mortality_cost_carbon ./ 1000) .* 1000  % we find 
    
    
end




% --------------------------------------------------------------------
% Values for UK report 2 - round to 100s
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valUKrep1 == 1


    % TOTAL emissions UK oil and gas fields (MtCO2e: add E6 to express as tCO2e)
    CO2_emissions_Jackdaw  = 19E6; 
    CO2_emissions_Rosebank = 122E6; 
 
    
    % Question 1: call function to compute the number of people facing one
    % extra heatwave due to the total emissions of the three fields
    clear valc_nr_children_facing_extra_heatwave_JackdawRosebank
    valc_nr_children_facing_extra_heatwave_JackdawRosebank(:,1)    = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_JackdawRosebank(:,2)    = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 5); % 5 is for heatwaves
    
    
    % include a small add-on table for the reference years 1960-1970
    clear valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref
    valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref      = valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref(end,:);                                                                                  % only store last row (decadal sum)
    valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref(2,:) = round(valc_nr_children_facing_extra_heatwave_JackdawRosebank(end,:) ./ valc_nr_children_facing_extra_heatwave_JackdawRosebank_ref(1,:) .* 100); % add a row with percentage change
        

    % Question 2: call function to compute the number of people facing one
    % extra ... due to the total emissions of the three fields
    % Droughts
    clear valc_nr_children_facing_extra_drought_JackdawRosebank
    valc_nr_children_facing_extra_drought_JackdawRosebank(:,1)     = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_JackdawRosebank(:,2)     = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 3); % 3 is for Droughts

    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_JackdawRosebank
    valc_nr_children_facing_extra_cropfailure_JackdawRosebank(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_JackdawRosebank(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 2); % 2 is for Crop Failure

    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_JackdawRosebank
    valc_nr_children_facing_extra_wildfire_JackdawRosebank(:,1)    = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_JackdawRosebank(:,2)    = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 1); % 1 is for wildfires

    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_JackdawRosebank
    valc_nr_children_facing_extra_tropcyclone_JackdawRosebank(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_JackdawRosebank(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 6); % 6 is for Tropical cyclones

    % River floods
    clear valc_nr_children_facing_extra_riverflood_JackdawRosebank
    valc_nr_children_facing_extra_riverflood_JackdawRosebank(:,1)  = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_JackdawRosebank(:,2)  = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 4); % 4 is for River floods

        
    % include a small add-on table for the reference years 1960-1970
    % Droughts
    clear valc_nr_children_facing_extra_drought_JackdawRosebank_ref
    valc_nr_children_facing_extra_drought_JackdawRosebank_ref(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_JackdawRosebank_ref(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_JackdawRosebank_ref      = valc_nr_children_facing_extra_drought_JackdawRosebank_ref(end,:);
    valc_nr_children_facing_extra_drought_JackdawRosebank_ref(2,:) = round(valc_nr_children_facing_extra_drought_JackdawRosebank(end,:) ./ valc_nr_children_facing_extra_drought_JackdawRosebank_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref
    valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref      = valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref(end,:);
    valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref(2,:) = round(valc_nr_children_facing_extra_cropfailure_JackdawRosebank(end,:) ./ valc_nr_children_facing_extra_cropfailure_JackdawRosebank_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref
    valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref      = valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref(end,:);
    valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref(2,:) = round(valc_nr_children_facing_extra_wildfire_JackdawRosebank(end,:) ./ valc_nr_children_facing_extra_wildfire_JackdawRosebank_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref
    valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref      = valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref(end,:);
    valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref(2,:) = round(valc_nr_children_facing_extra_tropcyclone_JackdawRosebank(end,:) ./ valc_nr_children_facing_extra_tropcyclone_JackdawRosebank_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % River floods
    clear valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref
    valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref(:,1) = mf_emissions2npeople(CO2_emissions_Jackdaw , TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref(:,2) = mf_emissions2npeople(CO2_emissions_Rosebank, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref      = valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref(end,:);
    valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref(2,:) = round(valc_nr_children_facing_extra_riverflood_JackdawRosebank(end,:) ./ valc_nr_children_facing_extra_riverflood_JackdawRosebank_ref(1,:) .* 100);
        
    
    
    % Question 3: compute the number of heat-related deaths between today and 2100 due to the total emissions of the three fields
    valc_mortality_JackdawRosebank(1,1) = floor(CO2_emissions_Jackdaw  ./ mortality_cost_carbon ./ 100) .* 100; % we find 
    valc_mortality_JackdawRosebank(1,2) = floor(CO2_emissions_Rosebank ./ mortality_cost_carbon ./ 100) .* 100 % we find 
    
    
end




% --------------------------------------------------------------------
% Values for UK report 2 - round to 100s
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valUKrep2 == 1


    % TOTAL emissions UK oil and gas fields (MtCO2e: add E6 to express as tCO2e)
    CO2_emissions_Penguins = 44E6; 

 
    
    % Question 1: call function to compute the number of people facing one
    % extra heatwave due to the total emissions of the three fields
    clear valc_nr_children_facing_extra_heatwave_Penguins
    valc_nr_children_facing_extra_heatwave_Penguins(:,1)    = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 5); % 5 is for heatwaves
    
    
    % include a small add-on table for the reference years 1960-1970
    clear valc_nr_children_facing_extra_heatwave_Penguins_ref
    valc_nr_children_facing_extra_heatwave_Penguins_ref(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_Penguins_ref      = valc_nr_children_facing_extra_heatwave_Penguins_ref(end,:);                                                                             % only store last row (decadal sum)
    valc_nr_children_facing_extra_heatwave_Penguins_ref(2,:) = round(valc_nr_children_facing_extra_heatwave_Penguins(end,:) ./ valc_nr_children_facing_extra_heatwave_Penguins_ref(1,:) .* 100);       % add a row with percentage change
        

    % Question 2: call function to compute the number of people facing one
    % extra ... due to the total emissions of the three fields
    % Droughts
    clear valc_nr_children_facing_extra_drought_Penguins
    valc_nr_children_facing_extra_drought_Penguins(:,1)     = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 3); % 3 is for Droughts

    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_Penguins
    valc_nr_children_facing_extra_cropfailure_Penguins(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 2); % 2 is for Crop Failure

    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_Penguins
    valc_nr_children_facing_extra_wildfire_Penguins(:,1)    = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 1); % 1 is for wildfires

    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_Penguins
    valc_nr_children_facing_extra_tropcyclone_Penguins(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 6); % 6 is for Tropical cyclones

    % River floods
    clear valc_nr_children_facing_extra_riverflood_Penguins
    valc_nr_children_facing_extra_riverflood_Penguins(:,1)  = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 2, 4); % 4 is for River floods

        
    % include a small add-on table for the reference years 1960-1970
    % Droughts
    clear valc_nr_children_facing_extra_drought_Penguins_ref
    valc_nr_children_facing_extra_drought_Penguins_ref(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_Penguins_ref      = valc_nr_children_facing_extra_drought_Penguins_ref(end,:);
    valc_nr_children_facing_extra_drought_Penguins_ref(2,:) = round(valc_nr_children_facing_extra_drought_Penguins(end,:) ./ valc_nr_children_facing_extra_drought_Penguins_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_Penguins_ref
    valc_nr_children_facing_extra_cropfailure_Penguins_ref(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_Penguins_ref      = valc_nr_children_facing_extra_cropfailure_Penguins_ref(end,:);
    valc_nr_children_facing_extra_cropfailure_Penguins_ref(2,:) = round(valc_nr_children_facing_extra_cropfailure_Penguins(end,:) ./ valc_nr_children_facing_extra_cropfailure_Penguins_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_Penguins_ref
    valc_nr_children_facing_extra_wildfire_Penguins_ref(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_Penguins_ref      = valc_nr_children_facing_extra_wildfire_Penguins_ref(end,:);
    valc_nr_children_facing_extra_wildfire_Penguins_ref(2,:) = round(valc_nr_children_facing_extra_wildfire_Penguins(end,:) ./ valc_nr_children_facing_extra_wildfire_Penguins_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_Penguins_ref
    valc_nr_children_facing_extra_tropcyclone_Penguins_ref(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_Penguins_ref      = valc_nr_children_facing_extra_tropcyclone_Penguins_ref(end,:);
    valc_nr_children_facing_extra_tropcyclone_Penguins_ref(2,:) = round(valc_nr_children_facing_extra_tropcyclone_Penguins(end,:) ./ valc_nr_children_facing_extra_tropcyclone_Penguins_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % River floods
    clear valc_nr_children_facing_extra_riverflood_Penguins_ref
    valc_nr_children_facing_extra_riverflood_Penguins_ref(:,1) = mf_emissions2npeople(CO2_emissions_Penguins, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 2, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_Penguins_ref      = valc_nr_children_facing_extra_riverflood_Penguins_ref(end,:);
    valc_nr_children_facing_extra_riverflood_Penguins_ref(2,:) = round(valc_nr_children_facing_extra_riverflood_Penguins(end,:) ./ valc_nr_children_facing_extra_riverflood_Penguins_ref(1,:) .* 100);
        
    
    
    % Question 3: compute the number of heat-related deaths between today and 2100 due to the total emissions of the three fields
    valc_mortality_Penguins(1,1) = floor(CO2_emissions_Penguins ./ mortality_cost_carbon ./ 100) .* 100; % we find 
    
    
end








% --------------------------------------------------------------------
% Values for Romania - round to 1000s
% only works if ms_valp is run !
% --------------------------------------------------------------------


if flags.valRomania == 1


    % TOTAL emissions UK oil and gas fields (MtCO2e: add E6 to express as tCO2e)
    CO2_emissions_NeptunDeep = 207E6; 

 
    
    % Question 1: call function to compute the number of people facing one
    % extra heatwave due to the total emissions of the three fields
    clear valc_nr_children_facing_extra_heatwave_NeptunDeep
    valc_nr_children_facing_extra_heatwave_NeptunDeep(:,1)    = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    
    
    % include a small add-on table for the reference years 1960-1970
    clear valc_nr_children_facing_extra_heatwave_NeptunDeep_ref
    valc_nr_children_facing_extra_heatwave_NeptunDeep_ref(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 5); % 5 is for heatwaves
    valc_nr_children_facing_extra_heatwave_NeptunDeep_ref      = valc_nr_children_facing_extra_heatwave_NeptunDeep_ref(end,:);                                                                             % only store last row (decadal sum)
    valc_nr_children_facing_extra_heatwave_NeptunDeep_ref(2,:) = round(valc_nr_children_facing_extra_heatwave_NeptunDeep(end,:) ./ valc_nr_children_facing_extra_heatwave_NeptunDeep_ref(1,:) .* 100);     % add a row with percentage change
        

    % Question 2: call function to compute the number of people facing one
    % extra ... due to the total emissions of the three fields
    % Droughts
    clear valc_nr_children_facing_extra_drought_NeptunDeep
    valc_nr_children_facing_extra_drought_NeptunDeep(:,1)     = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts

    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_NeptunDeep
    valc_nr_children_facing_extra_cropfailure_NeptunDeep(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure

    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_NeptunDeep
    valc_nr_children_facing_extra_wildfire_NeptunDeep(:,1)    = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires

    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_NeptunDeep
    valc_nr_children_facing_extra_tropcyclone_NeptunDeep(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones

    % River floods
    clear valc_nr_children_facing_extra_riverflood_NeptunDeep
    valc_nr_children_facing_extra_riverflood_NeptunDeep(:,1)  = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 2010, 2020, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods

        
    % include a small add-on table for the reference years 1960-1970
    % Droughts
    clear valc_nr_children_facing_extra_drought_NeptunDeep_ref
    valc_nr_children_facing_extra_drought_NeptunDeep_ref(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 3); % 3 is for Droughts
    valc_nr_children_facing_extra_drought_NeptunDeep_ref      = valc_nr_children_facing_extra_drought_NeptunDeep_ref(end,:);
    valc_nr_children_facing_extra_drought_NeptunDeep_ref(2,:) = round(valc_nr_children_facing_extra_drought_NeptunDeep(end,:) ./ valc_nr_children_facing_extra_drought_NeptunDeep_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Crop Failures
    clear valc_nr_children_facing_extra_cropfailure_NeptunDeep_ref
    valc_nr_children_facing_extra_cropfailure_NeptunDeep_ref(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 2); % 2 is for Crop Failure
    valc_nr_children_facing_extra_cropfailure_NeptunDeep_ref      = valc_nr_children_facing_extra_cropfailure_NeptunDeep_ref(end,:);
    valc_nr_children_facing_extra_cropfailure_NeptunDeep_ref(2,:) = round(valc_nr_children_facing_extra_cropfailure_NeptunDeep(end,:) ./ valc_nr_children_facing_extra_cropfailure_NeptunDeep_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Wildfires
    clear valc_nr_children_facing_extra_wildfire_NeptunDeep_ref
    valc_nr_children_facing_extra_wildfire_NeptunDeep_ref(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 1); % 1 is for wildfires
    valc_nr_children_facing_extra_wildfire_NeptunDeep_ref      = valc_nr_children_facing_extra_wildfire_NeptunDeep_ref(end,:);
    valc_nr_children_facing_extra_wildfire_NeptunDeep_ref(2,:) = round(valc_nr_children_facing_extra_wildfire_NeptunDeep(end,:) ./ valc_nr_children_facing_extra_wildfire_NeptunDeep_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % Tropical cyclones
    clear valc_nr_children_facing_extra_tropcyclone_NeptunDeep_ref
    valc_nr_children_facing_extra_tropcyclone_NeptunDeep_ref(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 6); % 6 is for Tropical cyclones
    valc_nr_children_facing_extra_tropcyclone_NeptunDeep_ref      = valc_nr_children_facing_extra_tropcyclone_NeptunDeep_ref(end,:);
    valc_nr_children_facing_extra_tropcyclone_NeptunDeep_ref(2,:) = round(valc_nr_children_facing_extra_tropcyclone_NeptunDeep(end,:) ./ valc_nr_children_facing_extra_tropcyclone_NeptunDeep_ref(1,:) .* 100);
        
    
    % include a small add-on table for the reference years 1960-1970
    % River floods
    clear valc_nr_children_facing_extra_riverflood_NeptunDeep_ref
    valc_nr_children_facing_extra_riverflood_NeptunDeep_ref(:,1) = mf_emissions2npeople(CO2_emissions_NeptunDeep, TCRE, exposure_perregion_BE, birth_years, 1960, 1970, GMT_BE, valp_cohort_size_abs, 1, 4); % 4 is for River floods
    valc_nr_children_facing_extra_riverflood_NeptunDeep_ref      = valc_nr_children_facing_extra_riverflood_NeptunDeep_ref(end,:);
    valc_nr_children_facing_extra_riverflood_NeptunDeep_ref(2,:) = round(valc_nr_children_facing_extra_riverflood_NeptunDeep(end,:) ./ valc_nr_children_facing_extra_riverflood_NeptunDeep_ref(1,:) .* 100);
        
    
    
    % Question 3: compute the number of heat-related deaths between today and 2100 due to the total emissions of the three fields
    valc_mortality_NeptunDeep(1,1) = floor(CO2_emissions_NeptunDeep ./ mortality_cost_carbon ./ 1000) .* 1000 
    
    
end






