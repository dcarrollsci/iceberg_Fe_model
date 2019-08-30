function ice = init_iceberg_size(L, dz)

%initialize iceberg size and shapes, based on length

%given L, outputs all other iceberg parameters
%dz : specify layer thickness desired, default is 10m

%Updated to make stable using Wagner et al. threshold, Sept 2017

%find keel depths, using barker or other method (auto, hotzel, mean)
K = keeldepth(L, 'barker'); % barker used in original paper

%default layer thickness, m
if nargin == 1
    
    dz = 10;
    
end

%iceberg coordinate system, x = length, y = width, z = thickness,

%now get underwater shape, based on Barker for K < 200, tabular for K > 200, and
%using L:W ratios of 1.62
ice = barker_carea(L, K, dz); %this gives you uwL, uwW, uwV, uwM, and vector Z down to keel depth

%from underwater volume, calculate above water volume
rhoi = 900; % (900 in original)
rati = rhoi ./ 1024; %ratio of ice density to water density
ice.totalV = (1/rati) .* nansum(ice.uwV,1)'; %~87% of ice underwater
ice.sailV = ice.totalV - nansum(ice.uwV,1)'; %sail volume is above water volume

%now figure out dimensions
ice.W = L./1.62; %iceberg width at waterline
ice.freeB = ice.sailV(:) ./ (L.*ice.W); %freeboard height
ice.L = L; %waterline length
ice.K = K;
ice.TH = ice.K + ice.freeB; %total thickness
ice.keeli = ceil(K ./ dz); %index of deepest iceberg layer
ice.dz = dz;
ice.dzk = -((ice.keeli-1) .* dz - K);

%now check if iceberg is stable
EC_thres = 0.92; %stability threshold from Wagner et al. 2017, if W/H < 0.92 then unstable
EC = ice.W ./ ice.TH;

%if less than threshold, then redo after fixing keel depth to make reasonable
if(EC < EC_thres)
    
    %method one, change keel depth to be shallower
    if(0)
        
        display(['Fixing keeldepth for L = ' num2str(L) ' m size class']);
        
        dbad = ice.TH - ice.W; %this is difference to get to stable thickness
        Knew = ice.K - rati .* dbad; %change by % of difference
        
        %redo calculations
        ice = barker_carea(L, Knew, dz); %this gives you uwL, uwW, uwV, uwM, and vector Z down to keel depth
        ice.totalV = (1/rati) .* nansum(ice.uwV,1)'; %~87% of ice underwater
        ice.sailV = ice.totalV - nansum(ice.uwV,1)'; %sail volume is above water volune
        ice.W = L ./ 1.62; %iceberg width at waterline
        ice.freeB = ice.sailV(:) ./ (L .* ice.W); % freeboard height
        ice.L = L; %waterline length
        ice.K = Knew; %keel thickness
        ice.TH = ice.K + ice.freeB; %total thickness
        ice.keeli = ceil(ice.K ./ dz); %index of deepest iceberg layer
        ice.dz = dz;
        ice.dzk = -((ice.keeli - 1) .* dz - ice.K);
        
        EC = ice.W ./ ice.TH;
        
    end
    
    %method two, change W to equal L, recalculate volumes
    if(1)
        
        display(['Fixing width to equal L, for L = ' num2str(L) ' m size class']);
        
        %use L:W ratio of to make stable, set so L:W makes EC=EC_thres
        Wtemp = EC_thres .* ice.TH;
        LWrat = floor(100 .* ice.L ./ Wtemp) ./ 100; %round down to hundreds place
        
        ice = barker_carea(L, K, dz, LWrat);
        
        %redo calculations
        ice.totalV = (1 ./ rati) .* nansum(ice.uwV,1)'; %~87% of ice underwater
        ice.sailV = ice.totalV - nansum(ice.uwV,1)'; %sail volume is above water volune
        ice.W = L ./ LWrat; %iceberg width at waterline
        ice.freeB = ice.sailV(:) ./ (L .* ice.W); %freeboard height
        ice.L = L; %waterline length
        ice.K = K; %keel thickness
        ice.TH = ice.K + ice.freeB; %total thickness
        ice.keeli = ceil(ice.K ./ dz); %index of deepest iceberg layer
        ice.dz = dz;
        ice.dzk = -((ice.keeli - 1) .* dz - ice.K);
        
        EC = ice.W ./ ice.TH;
        
    end
    
    if(EC < EC_thres)
        
        error('Still unstable, check W/H ratios')
        
    end
    
    
end