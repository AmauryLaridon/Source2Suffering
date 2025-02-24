  function A=hfbulktc(ur,zr,Ta,zt,rh,zq,Pa,Ts)
% HFBULKTC: computes sensible and latent heat fluxes and other variables.
% A=HFBULKTC(ur,zr,Ta,zt,rh,zq,Pa,Ts) computes the following variables:
%
%           Hs = sensible heat flux into ocean (W/m^2)
%           Hl = latent heat flux into ocean (W/m^2)
%           stress = wind stress (N/m^2)
%           usb = velocity friction scale (m/s)
%           Tsb = temperature scale (deg C)
%           Qsb = humidity scale (kg/kg)
%           L = Monin-Obukoff length (m)
%           zetu = zr/L
%           CD = drag coefficient
%           CT = temperature transfer coefficient (Stanton number)
%           CQ = moisture transfer coefficient (Dalton number)
%           RI = bulk Richardson number
%
% Based on the following buoy input data:
%
%           ur = wind speed (m/s) measured at height zr (m) 
%           Ta = air temperature (deg C) measured at height zt (m)
%           rh = relative humidity (%) measured at height zq (m)
%           Pa = air pressure (mb)
%           Ts = sea surface temperature (deg C)
%
% where ur, Ta, rh, Pa, and Ts may be either row or column vectors; zr,
% zt, and zq fixed scalars. rh and Pa may also be fixed scalars. Output
% variables are given as column vectors in A, i.e.,
%
%           A=[Hs Hl stress usb Tsb Qsb L zetu CD CT CQ RI]
%
% Code follows Edson and Fairall TOGA COARE code (version 2.0), modified 
% to include Rogers' weighting factor for unstable conditions.  Code does
% include gustiness, and assumes that the marine boundary layer height is
% known and constant over time for simiplicity. zr/L is limited to 
% be <=3.0 to ensure that the code converges to nonzero stress 
% and heat flux values for strongly stable 
% conditions.  The bulk Richardson number is computed between the sea 
% surface and zr as a diagnostic about whether turbulent boundary layer
% theory is applicable.  Code does not include either warm layer or cool
% skin effects to modify Ts.  See Fairall et al (1996), J. Geophys. Res.,
% 101, 3747-3764, for description of full TOGA COARE code and comparison
% with data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 28/Aug/98: Version 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP - 19/8/98 rewrote code quite a bit. Some inconsistencies in virtual
%      and real temperatures resolved (I hope). Loop structure rationalized. 
%      Computation of stress corrected to remove gusting. Gustiness
%      minimum setting done correctly.

% change column vectors to row vectors if necessary
[m n] = size(ur);
if m > n
  ur = ur';
end
[m n] = size(Ta);
if m > n
  Ta = Ta';
end
[m n] = size(rh);
if m > n
  rh = rh';
end
[m n] = size(Ts);
if m > n
  Ts = Ts';
end
[m n] = size(Pa);
if m > n
  Pa = Pa';
end

M=length(ur);

% create vectors for rh and Pa if scalars are input
if length(rh)==1 & M>1
  rh=rh*ones(1,M);
end
if length(Pa)==1 & M>1
  Pa=Pa*ones(1,M);
end


% Initialize various constants
as_consts;

tol=.001;    % Tolerance on Re changes to make sure soln has converged.

onethird=1./3;
o61=1/eps_air-1;   % 0.61 (moisture correction for temperature)



Le=(2.501-0.00237*Ts)*10^6;       % Latent heat
visc=viscair(Ta);                 % viscosity
Qsats=0.98*qsat(Ts);              % saturation specific humidity
                                  %    (reduced by 2% over salt water)
Q=(0.01.*rh).*qsat(Ta);           % specific humidity of air

T =Ta+CtoK;   % Convert to K
Ts=Ts+CtoK;

Tv=T.*(1 + o61*Q);                % Air virtual temperature
rho=(100*Pa)./(gas_const_R*Tv);   % air density


Dt=(T+0.0098.*zt)-Ts;     % adiabatic temperature difference
Dq=Q-Qsats;               % humidity difference

% compute initial neutral scaling coefficients
S=sqrt(ur.^2 + min_gustiness.^2);
cdnhf=sqrt(cdntc(S,zr,Ta)); % Smith's neutral cd as first guess

z0t=7.5*10^(-5);
ctnhf=kappa./log(zt./z0t);

z0q=z0t;
cqnhf=kappa./log(zq./z0q);

ust=cdnhf.*S;      % u_star_t (u_star including gustiness)
Ts =ctnhf.*Dt;     % T_star
Qs =cqnhf.*Dq;     % Q_star

% compute bulk Richardson number (as a diagnostic) - the "T"
% is probably not quite right - assumes T \approx Ts (good enough though).
RI=g.*zr.*(Dt + o61*T.*Dq)./(Tv.*S.^2);

Reu=0;Ret=0;Req=0;

% begin iteration loop to compute best ust, Ts, and Qs
for iter1=1:20;

    ReuO=Reu;RetO=Ret;ReqO=Req; % Save old values
    
    % Compute Monin-Obukov length (NB - definition given as eqn (7)
    % of Fairall et al (1996) probably wrong, following, e.g.
    % Godfrey and Bellars, JGR 96,22043-22048, 1991 and original code)
    bs=g*(Ts.*(1 + o61*Q) + o61*T.*Qs)./Tv; 
    L=(ust.^2)./(kappa*bs);
    % set upper limit on zr/L = 3.0 to force convergence under 
    % very stable conditions. Assume that zr, zt and zq comparable.
    L(L<zr/3 & L>0)=zr/3;
    
    zetu=zr./L;  % nondimensionalized heights
    zett=zt./L;
    zetq=zq./L;

    % surface roughness
    z0=(Charnock_alpha/g).*ust.^2 + R_roughness.*visc./ust;

    % Compute usb correction for non-neutral conditions
    cdnhf=kappa./(log(zr./z0)-psiutc(zetu));
    ust=cdnhf.*S;
  
    Reu=z0.*ust./visc;   % roughness Reynolds #
    [Ret,Req]=LKB(Reu);  % Compute other roughness reynolds #s

    % compute t and q roughness scales from roughness R#s
    z0t=visc.*Ret./ust;
    z0q=visc.*Req./ust;

    % compute new transfer coefficients at measurement heights
    cthf=kappa./(log(zt./z0t)-psittc(zett));
    cqhf=kappa./(log(zq./z0q)-psittc(zetq));

    % compute new values of Ts, Qs, and Tvsb
    Ts=cthf.*Dt;
    Qs=cqhf.*Dq;

    % estimate new gustiness
    Ws=ust.*(-CVB_depth./(kappa*L)).^onethird;
    wg=min_gustiness*ones(1,M); 
    j=find(zetu<0);                 % Convection in unstable conditions only
    wg(j)=max(min_gustiness,Beta.*Ws(j)); % Set minimum gustiness
    S=sqrt(ur.^2 + wg.^2);

%fprintf('%f %f %f\n',max(abs(Reu-ReuO)),max(abs(Ret-RetO)),max(abs(Req-ReqO)));

end % end of iteration loop
%plot([Reu' Ret' Req']);
%fprintf('%f %f %f\n',max(Reu),max(Ret),max(Req));

ii= abs(Reu-ReuO)>tol | abs(Ret-RetO)>tol | abs(Req-ReqO)>tol;
if any(ii),
 disp(['Algorithm did not converge for ' int2str(sum(ii)) ' values. Indices are:']);
 disp(find(ii)');
 warning('Not converged!');
end;

% compute fluxes into ocean
Hs=rho.*cp.*ust.*Ts;
Hl=rho.*Le.*ust.*Qs;

% compute transfer coefficients at measurement heights
CD=(ust./S).^2;
CT=ust.*Ts./(S.*Dt); % Stanton number
CQ=ust.*Qs./(S.*Dq); % Dalton number

% To compute mean stress we don't want to include the effects
% of gustiness which average out (in a vector sense).
stress=rho.*CD.*S.*ur;

% Output vector
A=[Hs' Hl' stress' ust' Ts' Qs'  L' zetu' CD' CT' CQ' RI'];





function y=psiutc(zet)
% PSIUTC: computes velocity profile function following TOGA/COARE.
% y=PSIUTC(zet) computes the turbulent velocity profile function given 
% zet = (z/L), L the Monin-Obukoff length scale, following Edson and
% Fairall TOGA COARE code (version 2.0) as modified to include Rogers' 
% weighting factor to combine the Dyer and free convection forms for 
% unstable conditions. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 28/Aug/98: version 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=length(zet);
c13=1.0./3.0;
sq3=sqrt(3.0);

% Stable conditions
y=-4.7*zet;

% Unstable conditions
j=find(zet<0);
zneg=zet(j);

% nearly stable (standard functions)
 x=(1-16.0.*zneg).^0.25;
 y1=2.0.*log((1+x)./2) + log((1+x.^2)./2) -2.*atan(x) + pi/2;

% Free Convective limit
 x=(1-12.87*zneg).^c13;
 y2=1.5*log((x.^2+x+1)./3) - sq3*atan((2.*x+1)/sq3) + pi/sq3;
	
% Weighted sum of the two
 F=1.0./(1+zneg.^2);
 y(j)=F.*y1+(1-F).*y2;




function y=psittc(zet)
% PSITTC: computes potential temperature profile following TOGA/COARE.
% y=PSITTC(zet) computes the turbulent potential temperature profile 
% function given zet = (z/L), L the Monin-Obukoff length scale, following 
% Edson and Fairall TOGA COARE code (version 2.0), as modified to use
% Rogers' weighting factor to combine the Dyer and free convective 
% forms for unstable conditions. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 28/Aug/98: version 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=length(zet);
c13=1.0./3.0;
sq3=sqrt(3.0);

% Stable conditions
y=-4.7.*zet;

% Unstable conditions
j=find(zet<0);
zneg=zet(j);

% Nearly stable (standard functions)
 x=(1-16.0.*zneg).^0.25;
 y1=2.0*log((1+x.^2)./2);

% Free Convective limit
 x=(1-12.87*zneg).^c13;
 y2=1.5.*log((x.^2+x+1)./3.0) - sq3.*atan((2.*x+1)./sq3) + pi./sq3;

% Weighted sum of the two
 F=1.0./(1+zneg.^2);
 y(j)=F.*y1 + (1-F).*y2;





function [Ret,Req]=LKB(Reu);
% LKB: computes Rougness reynolds numbers for temperature and humidity
%      using the theory of Liu, Katsaros and Businger, J. Atmos. Sci.
%      36, 1722-1735, 1979.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 28/Aug/98: version 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Ret=.177*ones(size(Reu));
   Req=.292*ones(size(Reu));
j=find(Reu>.11 & Reu<=.825);
   Ret(j)=1.376.*Reu(j).^0.929;
   Req(j)=1.808.*Reu(j).^0.826;
j=find(Reu>.825 & Reu<=3);
   Ret(j)=1.026./Reu(j).^0.599;
   Req(j)=1.393./Reu(j).^0.528;
j=find(Reu>3 & Reu<=10);
   Ret(j)=1.625./Reu(j).^1.018;
   Req(j)=1.956./Reu(j).^0.870;
j=find(Reu>10 & Reu<=30);
   Ret(j)=4.661./Reu(j).^1.475;
   Req(j)=4.994./Reu(j).^1.297;
j=find(Reu>30);
   Ret(j)=34.904./Reu(j).^2.067;
   Req(j)=30.790./Reu(j).^1.845;




