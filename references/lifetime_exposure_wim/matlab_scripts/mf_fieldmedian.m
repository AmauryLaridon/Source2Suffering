

% --------------------------------------------------------------------
% function to compute field averages
% works for multiple masks now
% --------------------------------------------------------------------


function [var_ap, varargout] = mf_fieldmedian(var, varargin)


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

        % get field median
        var_ap = nanmedian(var(:)   );
        var_mp = nanmedian(var(mask));       

    elseif length(size(var)) == 3

        % get field median
        nz     = size(var, 3);
        var_ap = NaN(nz, 1);
        var_mp = NaN(nz, 1);
        for j = 1:nz
            vari        = var(:,:,j);
            var_ap(j,1) = nanmedian(vari(:)   );
            var_mp(j,1) = nanmedian(vari(mask));
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

