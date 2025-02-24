

% --------------------------------------------------------------------
% function to plot data and the uncertainty bands
% --------------------------------------------------------------------


% note:
% - modified from http://stackoverflow.com/questions/17790817/displaying-trace-of-data-vs-time-with-band-of-uncertainity


function [h] = mf_plotuncertainty(x, y, stdev_up, stdev_down, color_line, color_range, lineprop, markersize, linewidth)


% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% initialise arrays
x = x(:);
y = y(:);


% % copy/past first and last column for smooth transitions at the figure boundaires
% x     = [x(1)-(x(2)-x(1)); x; x(end)+(x(end)-x(end-1))];
% y     = [y(end); y; y(1)];
% stdev = [stdev(end); stdev; stdev(1)];


% initalisation for transparency effect - doesn't work well yet!!
a = ones(size(x));



% --------------------------------------------------------------------
% visualisation
% --------------------------------------------------------------------


% % plot upper uncertainty range - variable transparency 
% patch('XData', [x; x(end:-1:1)], 'YData', [y + 1*stdev_up; y(end:-1:1)], ...
%       'FaceVertexAlphaData', [0*a; a], 'FaceAlpha', 'interp', 'EdgeColor', 'none','FaceColor', color_range); hold on;
% % plot lower uncertainty range - variable transparency                   
% patch('XData', [x; x(end:-1:1)], 'YData', [y - 1*stdev_down; y(end:-1:1)], ...
%       'FaceVertexAlphaData', [0*a; a], 'FaceAlpha', 'interp', 'EdgeColor', 'none','FaceColor', color_range); hold on;


% plot upper uncertainty range  
patch('XData', [x; x(end:-1:1)], 'YData', [y + 1*stdev_up; y(end:-1:1)], ...
      'FaceVertexAlphaData', [a; a], 'FaceAlpha', 'interp', 'EdgeColor', 'none','FaceColor', color_range); hold on;


% plot lower uncertainty range                   
patch('XData', [x; x(end:-1:1)], 'YData', [y - 1*stdev_down; y(end:-1:1)], ...
      'FaceVertexAlphaData', [a; a], 'FaceAlpha', 'interp', 'EdgeColor', 'none','FaceColor', color_range); hold on;

  
% now plot the data
h = plot(x, y, lineprop, 'markersize', markersize,'Color', color_line,'MarkerFaceColor', color_line, 'LineWidth', linewidth); hold on;


% plot the figure box
box on;


end

