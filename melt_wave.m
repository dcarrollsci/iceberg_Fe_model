function M = melt_wave(Wind_u, T_w, SeaIceC)
% Silva et al eqn for wave erosion term
% Mw = melt_wave(Wind_speed,T_surf, SeaIceConc)
%
% solves for melt rate Mw (in m/s), given
% T_w = surface temp of water
% Wind_u = in m/s (really relative to water speed, but assume Wind >>> water speeds)
% SeaIceC = sea ice concentration in % (0-1)
% 
% this is similar to mitberg wave formulation using the bigg option (after martin and adroft too),
% CIS model uses different formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SeaS = 1.5.*sqrt(Wind_u) + 0.1.*Wind_u; % sea state based on m/s winds
% BS = beaufort(Wind_speed); % get sea state given knots

IceTerm = 1 + cos(SeaIceC.^3 .* pi); % typo in their eqn. 7, should be + sign (See Gladstone et al. 2001)

M = (1/12).*SeaS.*IceTerm.*(T_w + 2); % in m/day

M = M ./ 86400; % now in m/s



%%%%% sub function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SeaState = beaufort(Wind_speed)
    % gives sea state given wind speed in knots
    
Wind_speed = abs(Wind_speed);
SeaState = 0;
if(Wind_speed > 1 && Wind_speed <= 3)
    SeaState = 1;
elseif(Wind_speed > 3 && Wind_speed <= 6)
    SeaState = 2;
elseif(Wind_speed > 6 && Wind_speed <= 10)
    SeaState = 3;
elseif(Wind_speed > 10 && Wind_speed <= 16)
    SeaState = 4;
elseif(Wind_speed > 16 && Wind_speed <= 21)
    SeaState = 5;
elseif(Wind_speed > 21 && Wind_speed <= 27)
    SeaState = 6;
elseif(Wind_speed > 27 && Wind_speed <= 33)
    SeaState = 7;
elseif(Wind_speed > 33 && Wind_speed <= 40)
    SeaState = 8;
elseif(Wind_speed > 40 && Wind_speed <= 47)
    SeaState = 9;
elseif(Wind_speed > 47 && Wind_speed <= 55)
    SeaState = 10;
elseif(Wind_speed > 55 && Wind_speed <= 63)
    SeaState = 11;
elseif(Wind_speed > 63)
    SeaState = 12;
end

if(length(Wind_speed)>1); SeaState = SeaState.*ones(length(Wind_speed),1); end
