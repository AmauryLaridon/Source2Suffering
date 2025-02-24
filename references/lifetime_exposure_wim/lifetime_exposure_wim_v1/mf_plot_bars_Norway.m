

% --------------------------------------------------------------------
% function to plot bar charts used in ETcHR expert opinion
% --------------------------------------------------------------------


function mf_plot_bars_Norway(exposure_ref, exposure_young, age_young, age_ref, extreme)




% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% define axes and sea color                               
darkcolor = [0.3  0.3  0.3 ]; % dark grey - 70% contrast (so 0.3) is advised
axcolor   = [0.5  0.5  0.5 ]; % normal grey
boxcolor  = [0.70 0.70 0.70]; % light grey
                           


% define colors
color_ref   = [102 104 177] ./ 256;
color_young = [177 102 103] ./ 256;





% --------------------------------------------------------------------
% manipulations
% --------------------------------------------------------------------


% get size of input field
[~, ncols] = size(exposure_ref);


if ncols == 3
    
    
    % include error bars
    exposure_bars_Norway =  [exposure_ref(1,2) exposure_young(1,2)];
    exposure_low_Norway  = ([exposure_ref(1,1) exposure_young(1,1)] - exposure_bars_Norway) .* -1;
    exposure_high_Norway =  [exposure_ref(1,3) exposure_young(1,3)] - exposure_bars_Norway;
    
    
elseif ncols == 1
    
    
    % do not include error bars
    exposure_bars_Norway =  [exposure_ref(1,1) exposure_young(1,1)];
   
    
end



% --------------------------------------------------------------------
% Visualisation
% --------------------------------------------------------------------


% create figure
% figure;
% set(gcf, 'color', 'w');
% set(gca,'color','w');
hold on;

% set bar width
barwidth = 0.4;
xlims    = [0.5 2.5];


% plot data
cats = categorical({age_ref, age_young});
cats = reordercats(cats,{age_ref, age_young}); 
ba = bar(cats, exposure_bars_Norway, barwidth, 'EdgeColor', 'none'); hold on;
ba.FaceColor = 'flat';
ba.CData(1,:) = color_ref;
ba.CData(2,:) = color_young;
set(gca, 'Fontsize', 13, 'Fontweight', 'Bold', 'Xcolor', boxcolor, 'Ycolor', boxcolor, 'box', 'off', 'TitleHorizontalAlignment', 'right');
xlabel('Birth years', 'Fontsize', 13, 'Fontweight', 'Bold'); 
ylabel('nr. people facing one extra extreme', 'Fontsize', 13, 'Fontweight', 'Bold'); 


% add error bar
if ncols == 3
er           = errorbar(cats, exposure_bars_Norway, exposure_low_Norway, exposure_high_Norway);    
er.Color     = darkcolor;                            
er.LineStyle = 'none';  
end


% % add exposure multiplication factors (EMF)
% text(ba(1).XEndPoints, ba(1).YEndPoints, {['\times' num2str(round(EMF_plot_OS(  ages==age_ref_plot)))],['\times' num2str(round(EMF_plot_OS(  ages==age_young)))]},'HorizontalAlignment','left','VerticalAlignment','middle', 'color', axcolor, 'Fontweight', 'Bold', 'fontsize', 14, 'rotation', 90)
% text(ba(2).XEndPoints, ba(2).YEndPoints, {['\times' num2str(round(EMF_plot_noOS(ages==age_ref_plot)))],['\times' num2str(round(EMF_plot_noOS(ages==age_young)))]},'HorizontalAlignment','left','VerticalAlignment','middle', 'color', axcolor, 'Fontweight', 'Bold', 'fontsize', 14, 'rotation', 90)


% Add title
title(extreme, 'color', axcolor, 'Fontweight', 'Bold', 'fontsize', 11)





end
