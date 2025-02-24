  function tau=stresslp(sp,z)
% STRESSLP: computes neutral wind stress following Large and Pond (1981).
% tau = STRESSLP(sp,z) computes the neutral wind stress given the wind
% speed at height z following Large and Pond (1981), J. Phys. Oceanog.,
% 11, 324-336. Assumes constant air density (1.22 kg/m^3). 
%
%   INPUT:   sp - wind speed (m/s)
%            z - measurement height (m)
%
%   OUTPUT:  tau - wind stress (N/m^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 26/8/98 - rewritten to use cdnlp

% Constant
as_consts;

[cd,u10]=cdnlp(sp,z);
tau=rho_air*(cd.*u10.^2);


