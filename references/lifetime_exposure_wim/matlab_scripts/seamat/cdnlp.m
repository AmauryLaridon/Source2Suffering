  function [cd,u10n]=cdnlp(sp,z)
% CDNLP: computes neutral drag coefficient following Large&Pond (1981).
% cd = CDNLP(sp,z) computes the neutral drag coefficient given the wind
% speed at height z following Large and Pond (1981), J. Phys. Oceanog.,
% 11, 324-336. Inputs can be scalars, vectors, or matrices.
%
%   INPUTS: sp - wind speed (m/s)
%           z - measurement height (m)
%
%   OUTPUT: cd_10 - neutral drag coefficient at 10m
%           u_10  - wind speed at 10m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorized RP 26/8/98

as_consts; % define physical constants


a=log(z./10)/kappa;  % Log-layer correction factor

tol=.001;       % Tolerance for iteration (m/s)


u10o=zeros(size(sp));
cd=1.15e-3*ones(size(sp));

u10n=sp./(1+a.*sqrt(cd));

% Notes - it appears to be slightly faster to do the math for all points
% then to extract submatrices of only those points that haven't yet
% converged in this iterative scheme

ii=abs(u10n-u10o)>tol;
while any(ii(:)),

  u10o=u10n;
  
  cd=(4.9e-4+6.5e-5*u10o);   % Compute cd_10(u_10)
  cd(u10o<10.15385)=1.15e-3;

  u10n=sp./(1+a.*sqrt(cd));   % Next iteration

  ii=abs(u10n-u10o)>tol; % Keep going until iteration converges
end;




 















