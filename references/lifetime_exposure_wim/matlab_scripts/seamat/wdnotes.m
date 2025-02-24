% WDNOTES: notes on estimating wave effects on wind measurements. 
% The following set of mfiles can be used to correct the wind speed 
% Ua measured at height za for the effects of the wave boundary layer
% following the empirical model presented by Large, Morzel, and 
% Crawford (1995), J. Phys. Oceanog., 25, 2959-2971. In particular,
% an analytic expression was found for the omega function (OMEGALMC)
% shown in their Fig. 9b, which allows the 'true' wind speed (Ut10) 
% and stress at 10m (assumed above the wave boundary layer height) 
% to be computed using WAVEDIST and the true wind speed (Uta) at the 
% measurement height za using WAVEDIS1. The Large et al model assumes 
% neutral stability (reasonable for high winds and wave conditions) 
% and uses a 10-m neutral drag law (CDNVE) based on Vera (1983;
% unpublished manuscript). This drag law follows Large and Pond (1982)
% for winds above 10 m/s but increases at lower wind speeds like 
% Smith (1987). The wave field is specified by the significant wave 
% height Hw.
%
% To compute 'true' wind speed Uta at za given Hw, use
%          Uta=wavedis1(Ua,za,Hw)   
% To compute 'true' wind speed Ut at 10m given Hw, use
%          [Ut10,(Ut10-U10)]=wavedist(Ua,za,Hw)
%
% Subroutines called:  
%          y=omegalmc(x)
%          cd10=cdnvera(u10)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

