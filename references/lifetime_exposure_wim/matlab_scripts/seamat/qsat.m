  function q=qsat(Ta,Pa)
% QSAT: computes specific humidity at saturation. 
% q=QSAT(Ta) computes the specific humidity (kg/kg) at satuation at
% air temperature Ta (deg C) using Tetens' formula for saturation vapor
% pressure from Buck (1981), J. App. Meteor., 1527-1532.  The dependence 
% on pressure is small (<0.5%) and has been removed using a mean pressure 
% of 1020 mb.  
%
%    INPUT:   Ta - air temperature  (deg C)
%             Pa - (optional) Pressure (mbars)
%
%    OUTPUT:  q  - saturation specific humidity  (kg/kg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin==1,
  as_consts;
  Pa=P_default; % pressure in mb
end;

a=(1.004.*6.112*0.6220)./Pa;
q=a.*exp((17.502.*Ta)./(240.97+Ta));


