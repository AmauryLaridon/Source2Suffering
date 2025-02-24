  function tau=stressve(sp,z)
% STRESSVE: computes stress using Vera (1983) neutral drag law.
% tau = STRESSVE(sp,z) computes the neutral wind stress given the wind
% speed at height z following Vera (1983) [see Large, Morzel, and Crawford
% (1995), J. Phys. Oceanog., 25, 2959-2971 (eqn. 8)]. Assumes z a fixed
% scalar, and constant air density (1.22 kg/m^3).
%
%   INPUT:  sp - wind speed (m/s)
%           z - measurement height (m)
%
%   OUTPUT: tau - wind stress (N/m^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 26/8/98 - rewritten to use cdnvera

% set constants
as_consts;

[cd,u10]=cdnve(sp,z);
tau=rho_air*(cd.*u10.^2);



