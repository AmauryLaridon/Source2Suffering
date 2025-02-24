

% --------------------------------------------------------------------
% function used to comute FLake profiles from default FLake output
% acknowledgement: Andrey Martynov, University of Bern
% --------------------------------------------------------------------


function [T_w_c H_ML_wm z_T_out] = mf_FLake_profiles(H_ML,T_wML,C_T,T_bot,nobs,depth_w_lk)



    % prepare for loop
    T_w_c   = nan(nanmax(depth_w_lk)+1,nobs);
    z_T_out = nan(nanmax(depth_w_lk)+1,nobs);
    
    % make function applicable to both spatial and temporal arrays
    if length(depth_w_lk) == 1
        depth_w_lk(1:nobs,1) = depth_w_lk;
    end
    
    % Calculate the water temperature profiles
    for i = 1:nobs
        
      % initialise evaluation site parameters
      z_T = [1E-3 1:depth_w_lk(i)]'; % Depth profile
      nz  = length(z_T);          % number of levels (e.g. 61 = 60+1)
        
      for j = 1:nz

        % fourth order polyniomial
        zeta = (z_T(j) - H_ML(i)) / (depth_w_lk(i) - H_ML(i));
        phi  = (40.0 / 3.0 * C_T(i) - 20.0 / 3.0) * zeta;
        phi  = phi + (18.0 - 30.0 * C_T(i)) * zeta.^2;
        phi  = phi + (20.0 * C_T(i) - 12.0) * zeta.^3;
        phi  = phi + (5.0 / 3.0 - 10.0 / 3.0 * C_T(i)) * zeta.^4;
        
        % assign temperatures
        if z_T(j) < H_ML(i)        % mixed layer
            T_w_c(j,i) = T_wML(i);
        else                       % thermocline
            T_w_c(j,i) = T_wML(i) - (T_wML(i)-T_bot(i)) * phi;
        end        
        
        % assign depth
        z_T_out(j,i) = z_T(j);
                
        % security
        if z_T(j) > depth_w_lk(i)
            T_w_c(j,i) = 0;
        end
        
      end
      
    end
    
    % get weekly running means of mixed layer depth
    windowSize = 24*7;
    H_ML_wm    = filter(ones(1,windowSize)/windowSize,1,H_ML);

end
