function M = melt_solar(Srad,alb)
% Melt from solar radiation in air, based on Condron's mitberg formulation, 
% would affect thickness above water only, assumes constant albedo for now
% 
% M = melt_solar(Srad,alb)
% 
% solves for melt rate M (in m/sec), given
% Srad: solar radiation flux downward (SW and LW), in W/m^2
% - note assumes iceberg albedo is 0.7
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define constants
Li = 3.33e5; % latent heat of fusion in ice J/kg
rho_i = 900; % density of ice
pera = 1 - alb; % percentage absorbed

% then melt rate in m/s
M = pera .* Srad ./ (rho_i .* Li);
