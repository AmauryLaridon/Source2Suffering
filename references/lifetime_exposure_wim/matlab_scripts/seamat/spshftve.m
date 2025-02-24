  function [sp2,ustar]=spshftve(sp1,z1,z2)
% SPSHFTVE: adjusts wind speed from z1 to z2 following Vera (1983).
% sp2 = SPSHFTVE(sp1,z1,z2) shifts the wind speed sp1 measured at z1 to
% z2 using the neutral drag law of Vera (1983) [see Large, Morzel, 
% and Crawford (1995), J. Phys. Oceanog., 25, 2959-2971 (eqn. 8)]. 
% Assumes z1 and z2 scalars.
%
%   INPUT:  sp1 - measured wind speed (m/s) 
%           z1 - measurement height (m)
%           z2 - desired height (m) of sp2
%
%   OUTPUT: sp2 - predicted wind speed (m/s)
%           ustar - friction velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RP 27/8/98 - completely rewritten to use cdnlp efficiently.

% set constants
as_consts;


% find cd and ustar
[cd10,sp10]=cdnve(sp1,z1);

ustar=sqrt(cd10).*sp10;

sp2=sp10+ustar.*log(z2./10)/kappa;

