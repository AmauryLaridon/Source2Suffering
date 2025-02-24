  function ew=satvap(Ta,P)
% SATVAP computes saturation vapor pressure
% q=satvap(Ta) computes the vapor pressure at satuation at
% air temperature Ta (deg C) (from Gill, Atmos-Ocean Dynamics, p606)
%
%    INPUT:   Ta- air temperature  (deg C)
%             p - pressure (mbars) (optional)
%
%    OUTPUT:  q  - saturation vapour pressure  (mbars)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 27/8/98 - small correction added (also from Gill)

as_consts;

if nargin<2,
  P=P_default;
end;
if nargin==0,
  Ta=Ta_default;
end;


ew=power(10,((0.7859+0.03477*Ta)./(1+0.00412*Ta)));

fw=1 + 1e-6*P.*(4.5+0.0006*Ta.^2);

ew=fw.*ew;
