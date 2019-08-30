function M = melt_buoyantwater(T_w, S_w, method)
% buoyant convection along sidewalls in water, based on bigg (condron)
% or CIS model formula
% 
% M = Melt_buoyantwater(T_w)
% (M = 7.62e-03 x dT + 1.3e-03 x dT^2)
%
% solves for melt rate M (in m/sec), given
% T_w = water temp. or dT if use CIS
% method = 'bigg' or 'cis', default is Bigg et al. model formula
% cis gives slightly larger melt rates, e.g., for 4 C ~0.01 m/day more
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
    method = 'bigg';
end

% define constants
Tf = -0.036 - (0.0499.*S_w) - (0.0001128.*S_w.^2); % freezing pt of seawater due to S changes 
Tfp = Tf.*exp(-0.19.*(T_w - Tf)); % fp temp

switch method
    case 'bigg'  % use bigg
        dT = T_w; 
        mday = 7.62e-3 .* dT + 1.3e-3.*dT.^2;
        
    case 'cis'  % use only barker
        dT = T_w - Tfp;
        mday = 7.62e-3 .* dT + 1.29e-3.*dT.^2;
end

% convert to m/s from m/day
M = mday ./ 86400; 