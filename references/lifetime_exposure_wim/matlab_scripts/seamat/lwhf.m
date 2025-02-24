  function qlw=lwhf(Ts,dlw,dsw)
% LWHF: computes net longwave heat flux following Dickey et al (1994).
% qlw=LWHF(Ts,dlw) computes the net longwave heat flux into the ocean.
%
% Following Dickey et al (1994), J. Atmos. Oceanic Tech., 11, 1057-1078,
% the incident longwave flux can be corrected for sensor heating due to
% insolation if you are using an Epply pygeometer. In this case, use
% qlw=LWHF(Ts,dlw,dsw).
%
% The surface emissivity is set at 0.98 based on Katsaros
% (1990), In: Surface Waves and Fluxes, ed. Geernaert and Plant, 672-701.
%
%   INPUT:  Ts  - sea surface temperature (deg C)
%           dlw - (measured) downward longwave flux (W/m^2)
%           dsw - (measured) insolation (W/m^2) (needed for Eppley pyrgeometers)
%
%   OUTPUT: qlw - net longwave heat flux (W/m^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Changes - RP 19/8/98 - remove correction for non-Epply pyrgeometers.

as_consts; % Get various constants

% convert degC to degK
ts=Ts+CtoK;

% correct dlw for insolation
% you only want to do this for Epply pyrgeometers!
if nargin==3,
  dlwc=dlw-0.036.*dsw;
else
  dlwc=dlw;
end;

% compute upward gray-body longwave flux
lwup=-emiss_lw.*sigmaSB.*(ts.^4);

% compute net flux
qlw=lwup + emiss_lw.*dlwc;



