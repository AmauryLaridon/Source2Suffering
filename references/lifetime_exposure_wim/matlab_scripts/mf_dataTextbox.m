
function obj = mf_dataTextbox(Xdata, Ydata, textstring, textcolor, boxcolor, BackgroundColor)

% based on
% https://nl.mathworks.com/matlabcentral/answers/346297-how-to-draw-an-arrow-using-non-normalized-coordinates


% call using e.g.
% mf_dataTextbox(Xdata,Ydata, agegroups{i,1}, axcolor, axcolor, seacolor)


% get axis handle
ax = gca;


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
% Ypixels = (Ydata - ax_ylim(1)) .* ax_pixels_per_ydata + axpos(2);   % convert data to normalised units
Ypixels = Ydata;                                                      % keep original (assume normalised)

length  = Xpixels(2) - Xpixels(1);
height  = Ypixels(2) - Ypixels(1);


% draw text box
obj = annotation('textbox',[Xpixels(1) Ypixels(1) length height], 'String', textstring, ...
      'EdgeColor', boxcolor, 'Color', textcolor, 'BackgroundColor', BackgroundColor,                           ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');


end