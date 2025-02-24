  function [cd,u10n]=cdnve(sp,z)
% CDNVE: computes neutral drag coefficient following Vera (1983).
% cd = CDNVE(sp10) computes the neutral drag coefficient given the wind
% speed at 10 m using the expression for friction velocity ustar derived 
% by E. Vera (1983) and published as eqn. 8 in Large, Morzel, and Crawford
% (1995), J. Phys. Oceanog., 25, 2959-2971. Range of fit to data is 1 to
% 25 m/s.
%
%   INPUT:  sp - wind speed (m/s)
%           z  - measurement height
%
%   OUTPUT: cd_10 - neutral drag coefficient at 10 m
%           u_10  - wind speed at 10 m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 26/8/98 - added reduction to 10m height.

% Constants in fit for drag coefficient
A=2.717e-3;
B=0.142e-3;
C=0.0764e-3;


as_consts;     % Other constant


a=log(z./10)/kappa;  % Log-layer correction factor
tol=.001;       % Tolerance for iteration (m/s)


u10o=zeros(size(sp))+.1;  % Don't start iteration at 0 to prevent blowups.
cd=(A./u10o + B + C*u10o);

u10n=sp./(1+a.*sqrt(cd));

% Notes - it appears to be slightly faster to do the math for all points
% then to extract submatrices of only those points that haven't yet
% converged in this iterative scheme

ii=abs(u10n-u10o)>tol;
while any(ii(:)),

  u10o=u10n;

  cd=(A./u10o + B + C*u10o);
  
  u10n=sp./(1+a.*sqrt(cd));   % Next iteration

  ii=abs(u10n-u10o)>tol; % Keep going until iteration converges
end;


