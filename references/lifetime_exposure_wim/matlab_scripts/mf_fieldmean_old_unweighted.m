

% --------------------------------------------------------------------
% function to compute field averages
% works for multiple masks now
% --------------------------------------------------------------------


function [var_ap, varargout] = mf_fieldmean(var, varargin)


% _ap: all pixels
% _mp: masked pixels


% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% Minimum number of land pixels that must contain a value in order to be
% retained
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
        var_ap = nanmean(nanmean(var));
        var_mp = nanmean(nanmean(var(mask)));

    elseif length(size(var)) == 3

        % get field mean for all pixel
        var_ap = nanmean(nanmean(var,1),2);

        % set array dimensions right
        var_ap = permute(var_ap,[3 2 1]);

        % get field mean for masked pixel
        var_mp = NaN(size(var_ap));
        for j = 1:length(var(1,1,:))
            vari        = var(:,:,j);
            var_mp(j,1) = nanmean(nanmean(vari(mask),1),2);
        end

    end
    
    
    % retain value only if data was available for more than X% of all land pixels in the subdomain 
    npixels   = length(find(       var(mask)));
    nobs      = length(find(~isnan(var(mask))));
    nobs_perc = (nobs / npixels) * 100;
    if nobs_perc < min_nobs;
        var_mp = NaN;
    end

    
    % store data
    varargout{i,:} = var_mp;

    
end

end

