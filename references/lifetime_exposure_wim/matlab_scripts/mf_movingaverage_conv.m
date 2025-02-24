

% --------------------------------------------------------------------
% function to apply a moving average filter
% input: vector x and window size w
% --------------------------------------------------------------------


% this function was inspired by this example
% http://matlabtricks.com/post-11/moving-average-by-convolution



function y = mf_movingaverage_conv(x, w)

   k = ones(1, w) / w;
   y = conv(x, k, 'same');

end

