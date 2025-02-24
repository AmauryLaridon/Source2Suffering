
function obj = mf_dataArrow(Xdata,Ydata, ax, color, LineWidth, HeadStyle)


% based on
% https://nl.mathworks.com/matlabcentral/answers/346297-how-to-draw-an-arrow-using-non-normalized-coordinates


% get axis handle
% ax = gca;


%get axes drawing area position in pixels relative to figure
oldunits = get(ax, 'Units');
set(ax, 'Units', 'Normalized');
axpos = get(ax, 'Position');
set(ax, 'Units', oldunits);


%get axes drawing area in data units
ax_xlim = xlim(ax);
ax_ylim = ylim(ax);
ax_pixels_per_xdata = axpos(3) ./ diff(ax_xlim);
ax_pixels_per_ydata = axpos(4) ./ diff(ax_ylim);


%these are figure-relative
Xpixels = (Xdata - ax_xlim(1)) .* ax_pixels_per_xdata + axpos(1);
Ypixels = (Ydata - ax_ylim(1)) .* ax_pixels_per_ydata + axpos(2);


% draw arrow
obj = annotation('arrow', Xpixels, Ypixels, 'Units', 'pixels', 'color', color, 'LineWidth', LineWidth, 'HeadStyle', HeadStyle);



end