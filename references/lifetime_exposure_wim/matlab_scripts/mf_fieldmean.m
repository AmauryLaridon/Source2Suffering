

% --------------------------------------------------------------------
% function to compute field averages weighted by pixel area
% works for multiple masks now
% --------------------------------------------------------------------


function [var_ap, varargout] = mf_fieldmean(var, area, varargin)


% _ap: all pixels
% _mp: masked pixels


% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% Minimum number of land pixels that must contain a value to be retained
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
        var_ap = nansum(var       .* area      , [1 2]) ./ nansum(area       , [1 2]);
        var_mp = nansum(var(mask) .* area(mask), [1 2]) ./ nansum(area(mask) , [1 2]);

        
    elseif length(size(var)) == 3
        
        
        % get number of steps in third dimension
        nz = size(var,3);

        
        if    length(size(area)) == 2 
        
            
            % get field mean
            var_ap = squeeze( nansum(var .* repmat(area        , 1, 1, nz) , [1 2] ) ./ repmat( nansum( area      , [1 2]), 1, 1, nz) );
            var_mp = squeeze( nansum(var .* repmat(area .* mask, 1, 1, nz) , [1 2] ) ./ repmat( nansum( area(mask), [1 2]), 1, 1, nz) );

            
        elseif length(size(area)) == 3
            
            % get field mean
            var_ap = squeeze( nansum(var .* area                          , [1 2] ) ./ nansum( area                          , [1 2]));
            var_mp = squeeze( nansum(var .* area .* repmat(mask, 1, 1, nz), [1 2] ) ./ nansum( area .* repmat(mask, 1, 1, nz), [1 2]));
                        
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

