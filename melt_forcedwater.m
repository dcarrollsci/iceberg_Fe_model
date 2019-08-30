function [M, T_sh, T_fp] = melt_forcedwater(T_far, S_far, P_base, U_rel)
% Silva et al eqn, using parameters from Holland and Jenkins
% M = melt_forcedwater(T_far,S_far,P_base,U_rel)
%
% solves for melt rate M (in m/sec), given
% T_far = farfield temp
% S_far = farfield S
% P_base = pressure at base
% U_rel = water speed relative to iceberg surface (this could be the ambient velocity
%             reported by Jenkins' plume model, or a horizontal relative velocity moving past iceberg
%
% also reports back Tsh and T_fp which is difference between T_far and local fp (Tsh = T_far - (aS_far + b + cP_base))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants (these all same as in Jenkins now, for comparison)
a = -5.73e-2; % S contribution to freezing temp 
b = 8.32e-2;  % constant
c = -7.61e-4; % pressure contribution in C/dbar
% so T_fp = aS_base + b + cP_base (where S_base,P_base are at base of iceberg)

% changing these Gamma's to Jenkins made fit better, especially at depth
GT = 1.1e-3; %6e-4; %second value in Silva, first in Jenkins model  % heat transfer coefficient
GS = 3.1e-5; %2.2e-5; %second value in Silva, first in Jenkins model % salt transfer coefficient
L = 3.35e5;  % latent heat of fusion of ice (in J/kg)
cw = 3974;   % specific heat of water (J/kg/C)
ci = 2009;   % " of ice
DT = 15;     % temp difference between iceberg core and bottom surface (this is for Antarctica, maybe less for Greenland?)

% solve eqns 3,4,5 for M and you get a quadratic with lots of crazy terms
%    that depend on farfield temp, salinity, and pressure at base, as well 
%    as abs(u) which is velocity magnitude of relative motion of iceberg to water

T_fp = a.*S_far + b + c.*P_base; % fp temp
T_sh = T_far - T_fp; % temp above fp 
% terms in quadratic
A = (L + DT.*ci)./(U_rel.*GT.*cw);  % term on quadratic
B = -((L + DT.*ci).*GS./(GT.*cw) - a.*S_far - T_sh);
C = -U_rel.*GS.*T_sh; 

ROOT = B.^2 - 4.*A.*C; 
ROOT(ROOT<0) = nan; 

% find roots of quadratic
Mtemp1 = (-B + sqrt(ROOT)) ./ (2.*A);
Mtemp2 = (-B - sqrt(ROOT)) ./ (2.*A);

M = min(Mtemp1,Mtemp2); % this is in m/sec
% now get rid of melt rates below fp, i.e., where T_sh<0
bad = T_sh < 0; % will be zero where no thermal driving happening
M(bad) = 0; 
M(isnan(M)) = 0;

% make melt rates positive
M = -M; 