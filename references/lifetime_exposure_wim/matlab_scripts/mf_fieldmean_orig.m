

% --------------------------------------------------------------------
% function to compute field averages weighted by pixel area
% works for multiple masks now
% --------------------------------------------------------------------


function [var_ap, varargout] = mf_fieldmean_orig(var, area, varargin)


% _ap: all pixels
% _mp: masked pixels


% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% Minimum number of land pixels that must contain a value in order to be retained
min_nobs = 50; % [%]



% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------


% loop over masks
for i=1:length(varargin)

    
    % set the mask
    mask = varargin{i};


    if    length(size(var)) == 2

        % get field mean
        var_ap = nansum(nansum(var       .* area      )) ./ nansum(nansum(area      ));
        var_mp = nansum(nansum(var(mask) .* area(mask))) ./ nansum(nansum(area(mask)));

    elseif length(size(var)) == 3

        % get field mean for all and masked pixels
        var_ap = NaN(size(var,3),1);
        var_mp = NaN(size(var,3),1);
        for j = 1:size(var,3)
            vari        = var(:,:,j);
            var_ap(j,1) = nansum(nansum(vari       .* area      )) ./ nansum(nansum(area      ));
            var_mp(j,1) = nansum(nansum(vari(mask) .* area(mask))) ./ nansum(nansum(area(mask)));
        end

    end
    
    
    % retain value only if data was available for more than X% of all land pixels in the subdomain 
    npixels   = length(find(       var(mask)));
    nobs      = length(find(~isnan(var(mask))));
    nobs_perc = (nobs / npixels) * 100;
    if nobs_perc < min_nobs
        var_mp = NaN;
    end

    
    % store data
    varargout{i,:} = var_mp;

    
end

end

