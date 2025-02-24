

% --------------------------------------------------------------------
% function to bin data 
% --------------------------------------------------------------------


function [binmean, binmedian, binstd, q25, q75, q90, q99, q999] = mf_bin(time_series, bin, nbins)


% if bin is an empty variable, then data is binned according to itself


% assign bins if necessary
if isempty(bin)   
    
    % get bins (use "tabulate(bin)" for some bin information)
    bin = ceil(nbins * tiedrank(time_series) / length(time_series));
    
end


% bin variables according to 'bin'
binmean   = NaN(nbins,1);
binmedian = NaN(nbins,1);
binstd    = NaN(nbins,1);
q25       = NaN(nbins,1);
q75       = NaN(nbins,1);
q90       = NaN(nbins,1);
q99       = NaN(nbins,1);
q999      = NaN(nbins,1);
for i=1:nbins
    binmembers   = time_series(bin == i);
    binmean(i)   = nanmean(binmembers);
    binmedian(i) = nanmedian(binmembers);
    binstd(i)    = nanstd (binmembers);
    q25(i)       = quantile(binmembers,0.25);
    q75(i)       = quantile(binmembers,0.75);
    q90(i)       = quantile(binmembers,0.90);
    q99(i)       = quantile(binmembers,0.99);
    q999(i)      = quantile(binmembers,0.999);
end
binmean(isnan(binmean)) = 0;
binstd (isnan(binstd )) = 0;


end

