

% --------------------------------------------------------------------
% function to plot regional averages under current pledges
% --------------------------------------------------------------------


function mf_plot_pie(EMF_perregion_NDC, EMF_perregion_20, EMF_perregion_15, ages, ages_pie, ind_age, regions, extremes, ind_extreme, ind_region, panel, color_NDC, color_20, color_15, EMF_label_step, EMF_label_max, flag_labels)




% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% define axes and sea color                               
axcolor  = [0.5  0.5  0.5 ]; % normal grey
                           

% define the origin
origin = [0, 0];


% define rotation of EMF labels (don't take pi/4 because then aligns with South Asia)
EMF_label_rotation = 3*pi/8;


% define angles of the pieces of the pie
cohort_size_rad = cell2mat(regions.cohort_size_rel(ind_region,ind_age)) .* 2.*pi; % width of each piece of the pie
cohort_size_rad = [0; cumsum(cohort_size_rad)] + pi/2;                            % convert to cumulative for plotting adjacent on a circle and rotate by 90° countreclockwise to have 


% get number of extremes
nextremes = size(EMF_perregion_NDC,1) - 1;



% --------------------------------------------------------------------
% Visualisation
% --------------------------------------------------------------------


% create figure
% figure;
% set(gcf, 'color', 'w');
% set(gca,'color','w');
hold on;


% plot pieces of the pie
for i=1:length(ind_region)
    P1 = mf_plot_arc(cohort_size_rad(i), cohort_size_rad(i+1), origin, EMF_perregion_NDC(ind_extreme, ind_region(i), ages==ages_pie(ind_age)), color_NDC);
    P2 = mf_plot_arc(cohort_size_rad(i), cohort_size_rad(i+1), origin, EMF_perregion_20( ind_extreme, ind_region(i), ages==ages_pie(ind_age)), color_20 );
    P3 = mf_plot_arc(cohort_size_rad(i), cohort_size_rad(i+1), origin, EMF_perregion_15( ind_extreme, ind_region(i), ages==ages_pie(ind_age)), color_15 );
    set(P1,'edgecolor', color_NDC,'linewidth', 2); % draw contour around biggest piece of the pie
end


% plot concentric circles with EMF as radius
for radius=EMF_label_step:EMF_label_step:EMF_label_max     % EMF radius    

    % plot concentric circles with EMF as radius
    rectangle('Position', [origin-radius 2*radius 2*radius],'Curvature',[1 1], 'Edgecolor', axcolor)

% end
% for radius=2*EMF_label_step:2*EMF_label_step:EMF_label_max     % EMF radius    
    
    % add EMF label
    text(origin(1) + radius*cos(EMF_label_rotation), origin(2) + radius*sin(EMF_label_rotation), ['\times' num2str(radius)], 'color', axcolor , 'Fontsize', 9, 'Fontweight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')


end
axis equal
axis off


% plot region labels
if flag_labels == 1
    pieHandle = pie(cell2mat(regions.cohort_size_rel(ind_region,ind_age)), regions.name_2lines(ind_region));
    pText = findobj(pieHandle,'Type','text');
    for iHandle = 2:2:2*numel(regions.name(ind_region))
        pieHandle(iHandle).Position = EMF_label_max * pieHandle(iHandle).Position;
        pieHandle(iHandle).Color    = axcolor;
    end
    for iHandle = 1:2:2*numel(regions.name(ind_region))
        pieHandle(iHandle).FaceAlpha = 0;
        pieHandle(iHandle).EdgeAlpha = 0;
    end
end


% get axes limits
ylims = ylim;
xlims = xlim;


% % add panel letter as text
% text(xlims(1) - EMF_label_step, ylims(2) , [panel ' ' num2str(ages_pie(ind_age)) ' yr'],'ver','bottom','hor','left','Fontsize', 15, 'Fontweight', 'bold', 'color', axcolor)


% % load and plot pictograms
% if ind_extreme == nextremes+1
%     pictogram = imadjust(imcomplement(imread('pictogram_all_2x3.png')),[0.3 0.3 0.3; 1 1 1],[]);
%     ax2 = axes('Position',[0.75 0.82 0.25 0.17]); %#ok<NASGU>
% else
%     pictogram = imadjust(imcomplement(imread(['pictogram_' extremes{ind_extreme} '.png'])),[0.3 0.3 0.3; 1 1 1],[]);
%     ax2 = axes('Position',[0.90 0.90 0.10 0.10]); %#ok<NASGU>
% end
% imagesc(pictogram); hold off   % add pictogram
% axis('off', 'image');




end
