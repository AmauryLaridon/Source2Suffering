

% --------------------------------------------------------------------
% function to apply a 1D moving average filter to 2D data, 
% and assign it to the last year of each window 
% input: vector x and window size w
% --------------------------------------------------------------------


% the value of the year 2010 represents the average during 2001-2010; cfr. Seneviratne et al., Nature 2016)



function out = mf_movingaverage(in, w, in_init)



% paste earlier values to matrix to avoid getting NaN at the start of you calculation
if ~isempty(in_init)
    in = [in_init(end-w+2:end,:); in];
end


% prepare for loop
[ny, nx] = size(in);
out      = NaN(ny, nx);


% loop over points in x-direction
for j=1:nx; 
    
    % loop over points in y-direction
    for i=w:ny;
        
        % get moving average in y-direction
        out(i,j) = nanmean(in(i-w+1:i, j));
        
    end
    
end


% remove these values again
if ~isempty(in_init)
    out = out(w:end,:);
end


end

