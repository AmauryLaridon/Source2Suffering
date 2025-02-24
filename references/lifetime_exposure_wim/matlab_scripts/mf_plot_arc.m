function P = mf_plot_arc(a, b, origin, r, face_color)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians, 
% b is end of arc in radians, 
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4, [9, -4],3)
% Author:  Matt Fig
% source: https://nl.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle
h = origin(1);
k = origin(2);

t = linspace(a,b);
x = r*cos(t) + h;
y = r*sin(t) + k;
x = [x h x(1)];
y = [y k y(1)];
P = fill(x,y,face_color);
set(P,'edgecolor','none')
% axis([h-r-1 h+r+1 k-r-1 k+r+1]) 
% axis square;
if ~nargout
    clear P
end