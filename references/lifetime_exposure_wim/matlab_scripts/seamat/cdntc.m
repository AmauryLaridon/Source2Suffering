  function [cd,u10]=cdntc(sp,z,Ta)
% CTDTC: computes the neutral drag coefficient following Smith (1988).
% cd = CDNTC(sp,z,Ta) computes the neutral drag coefficient given the 
% wind speed and air temperature at height z following Smith (1988),
% J. Geophys. Res., 93, 311-326. Assumes sp and Ta are both column 
% or row vectors and z a fixed scalar.
%
%   INPUT:  sp - wind speed (m/s)
%           z - measurement height (m)
%           Ta - air temperature (deg C) (optional)
%
%   OUTPUT: cd_10 - neutral drag coefficient at 10m
%           u_10  - wind speed at 10m (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorized RP 26/8/98

as_consts; % Define constants

if nargin==2,
  Ta=Ta_default;
end;

% Iteration endpoint
tol=.00001;


visc=viscair(Ta);

% remove any sp==0 to prevent division by zero
i=find(sp==0);
sp(i)=.1.*ones(length(i),1);


% initial guess
ustaro=zeros(size(sp));
ustarn=.036.*sp;

% iterate to find z0 and ustar

ii=abs(ustarn-ustaro)>tol;
while any(ii(:)),

  ustaro=ustarn;
  z0=Charnock_alpha.*ustaro.^2./g + R_roughness*visc./ustaro;
  
  ustarn=sp.*(kappa./log(z./z0));
 
  ii=abs(ustarn-ustaro)>tol;
end

sqrcd=kappa./log((10)./z0);
cd=sqrcd.^2;

u10=ustarn./sqrcd;



