% beautiful matlab boxplots
% John Bogovic
% November 2011
% http://www.iacl.ece.jhu.edu/JohnMatlabTips

%% some parameters for the plots

% the order that the tags should appear on the boxplots
grpOrder = {'Group1A','Group1B', 'Group1C', ...
    'Group2A','Group2B', 'Group2C'};

% positions on the x-axis to drop boxplots
pos = [1 2 3 7 8 9];  

% this tag tells matlab which boxes should be the same color
% let's color the lettered groups the same way.
allcolorgroups = {'A','B','C','A','B','C'};

% a color map
cmap = hsv2rgb([0 0.6 0.6; 0.3 0.6 0.6; 0.6 0.6 0.6]);

%% gen some fake data

datmtx = randn(600,1);

% randomly pick a group for each measurement
grpi = randi(6,600,1);
grp={grpOrder{grpi}};
% make an appropriate colorgrouping
cgrp = {allcolorgroups{grpi}};

%% make boxplots
figure;

% make the boxplots
% the 'compact' plotstyle goes a long way toward making things nice
% group colors, and spacings along the x-axis are great in guiding the eye
boxplot(datmtx,grp,'plotstyle','compact','grouporder',grpOrder,...
    'symbol','k.','positions',pos,'colorgroup',cgrp,'colors',cmap,...
    'medianstyle','line');

grid on;

% make median lines black and big
set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',2);

% make outlier dots gray and big
set(findobj(gcf,'Tag','Outliers'),'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerSize',5);