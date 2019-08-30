function M = melt_forcedair(T_air, U_rel, L)
% Forced convection in air, based on Condron's mitberg formulation, 
% this is really M = HF / (rho_i x L_i), where HF is heat flux
% 
% M = Melt_forcedair(T_air, U_rel, L)
% 
% solves for melt rate M (in m/sec), given
% T_air = temp in air (assumes temp in ice = -4)
% U_rel = air speed relative to iceberg (typically just wind speed, as these are >> iceberg speeds)
% L: iceberg L (m)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define constants
Tice = -4; 
Li = 3.33e5; % latent heat of fusion in ice J/kg
rho_i = 900; % density of ice
visca = 1.46e-5; %kinematic viscosity of air
diffa = 2.16e-5; %thermal diffusivity of air
conda = 0.0249; % thermal conductivty of air
cold = T_air<0; % make these M's = 0

% dimensionless #'s
Pr = visca ./ diffa; % Prandtl number
Re = abs(U_rel).*L ./ visca; %Reynolds number based on relative air speed
Nu = 0.058 .* (Re.^0.8) ./ (Pr.^0.4); %Nusselt number

% then heat flux
HF = (1./L) .* (Nu .* conda .* (T_air - Tice));

% then melt rate in m/s
M = HF ./ (rho_i .* Li);

M(cold) = 0; % any T_air's < 0 make melt rate = 0; 