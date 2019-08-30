function M = melt_wave_bigg(Wind_u)
% Bigg et al eqn for wave erosion term
% Mw = melt_wave(Wind_speed)
%
% here, melt rate is only function of wind speed, no SeaIce or T dependence
% Wind_u = in m/s (really relative to water speed, but assume Wind >>> water speeds)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to knots
Wind_u = 1.94 .* Wind_u; % knots are 2x m/s

BS = beaufort(Wind_u); % get sea state given knots

M = (1/2) .* BS; % in m/day

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
