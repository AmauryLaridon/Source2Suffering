
% --------------------------------------------------------------------
% function to draw a circle in a matrix using the integer midpoint 
% circle algorithm. Does not miss or repeat pixels.
% Created by : Peter Bone
% Created : 19th March 2007
% --------------------------------------------------------------------


function i = mf_MidpointCircle(i, radius, xc, yc, value, flag_filled)


% Make sure you are working with integers
xc = int16(xc);
yc = int16(yc);

x = int16(0);
y = int16(radius);
d = int16(1 - radius);


% assign value to pixels on x- and y-axes
i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;


% start while loop for other points
while ( x < y - 1 )
    x = x + 1;             % move up along positive x-axis
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    i( x+xc,  y+yc) = value; % for the first octant
    i( y+xc,  x+yc) = value; % for the other octants... (mirrored)
    i( y+xc, -x+yc) = value;
    i( x+xc, -y+yc) = value;
    i(-x+xc, -y+yc) = value;
    i(-y+xc, -x+yc) = value;
    i(-y+xc,  x+yc) = value;
    i(-x+xc,  y+yc) = value;
end



if flag_filled == 1
    
    % fill the circle
    for ii = xc-int16(radius):xc+(int16(radius))
        for jj = yc-int16(radius):yc+(int16(radius))
            tempR = sqrt((double(ii) - double(xc)).^2 + (double(jj) - double(yc)).^2);
            if(tempR <= double(int16(radius)))
                i(ii,jj)=value;
            end
        end
    end

end


end

