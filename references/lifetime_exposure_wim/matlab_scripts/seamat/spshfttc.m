  function [sp2,ustar]=spshfttc(sp1,z1,z2,Ta)
% SPSHFTTC: adjusts wind speed from z1 to z2 following Smith (1988).
% sp2 = SPSHFTTC(sp1,z1,z2,Ta) shifts the wind speed sp1 measured at z1 to
% z2 using the neutral drag coefficient given the wind speed and air
% temperature at height z following Smith (1988), J. Geophys. Res., 93,
% 311-326. Assumes z1 and z2 scalars. Ta may be a constant. 
%
%   INPUT:  sp1 - measured wind speed (m/s)
%           z1 - measurement height (m)
%           z2 - desired height (m)
%           Ta - air temperature (deg C) (optional)
%
%   OUTPUT: sp2 - predicted wind speed 
%           ustar - fiction velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 27/8/98 - completely rewritten to use cdntc efficiently.

% set constants
as_consts;

if nargin==3,
  Ta=Ta_default;
end;

% find cd and ustar
[cd,sp10]=cdntc(sp1,z1,Ta);

ustar=sqrt(cd).*sp10;

sp2=sp10+ustar.*log(z2./10)/kappa;
