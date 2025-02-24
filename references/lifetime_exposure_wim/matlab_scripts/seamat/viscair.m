  function  v=viscair(Ta)
% VISCAIR: computes viscosity of air 
% v=VISCAIR(Ta) computes the kinematic viscosity of dry air (m^2/sec)
% as a function of Ta (deg C) following Andreas (1989), CRREL Report 
% 89-11.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=(1+(6.542.*Ta.*10^(-3))+(8.301.*Ta.^2.*10^(-6))-(4.84.*Ta.^3.*10^(-9)));
v=(1.326.*10^(-5)).*v;

