  function tau=stresstc(sp,z,Ta)
% STRESSTC: computes the neutral wind stress following Smith (1988).
% tau = STRESSTC(sp,z,Ta) computes the neutral wind stress given the 
% wind speed and air temperature at height z following Smith (1988),
% J. Geophys. Res., 93, 311-326. Assumes a constant air density (1.22 kg/m^3).
%
%   INPUT:  sp - wind speed (m/s)
%           z  - measurement height (m)
%           Ta - air temperature (deg C) (optional)
%
%   OUTPUT: tau - wind stress (N/m^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 26/8/98 Rewritten to use cdntc

% compute constants
as_consts;

if nargin==2,
  Ta=Ta_default;
end;

[cd,u10]=cdntc(sp,z,Ta);

tau=rho_air*(cd.*u10.^2);

