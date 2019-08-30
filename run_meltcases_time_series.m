clear
close all;

dataDir = '/Users/carrolld/Documents/research/greenland/iceberg_iron/raw_data/all/';

A = xlsread([dataDir 'Iceberg data compilation complete.xlsx']);
speckleFe = A(1:216,6) .* 10^-9;
speckleFe(isnan(speckleFe)) = [];

saveDir = '/Users/carrolld/Documents/research/greenland/iceberg_iron/mat/moon_iceberg_model/';

addpath ./Icebergs_code

%%

load('./Infiles/Seasonal/TS_monthly.mat');
load('./Infiles/Seasonal/monthly_MET.mat');

%% 

%number of icebergs
%L = 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000
per = [0.47; 0.39; 0.08; 0.03; 0.01; 0.01; 0.004; 0.002; 0.001; 7e-4; 5e-4; ...
    3e-4; 2e-4; 2e-4; 1e-4; 1e-4; 5e-5; 5e-5; 5e-5; 5e-5];

N = 24868; % from Dan Sulak, total # of icebergs (median value)
Nbergs = N .* per(:); %number of icebergs

%L = 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000
sice = [0.07; 0.03; 0.03; 0.03; 0.03; 0.01; 0.01; 0.01; 0.01; 0.01; 0.005; 0.005; ...
    0.005; 0.005; 0.005; 0.005; 0.005; 0.001; 0.001; 0.0005];

Nsice = sice .* Nbergs;
Nbergs_adj = Nbergs - Nsice;

%monthly sea ice concentration. (or total ice cover)
%SEP, OCT, NOV, DEC, JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG
MONcice = [0.35 0.3 0.36 0.5 0.65 0.7 0.7 0.6 0.55 0.45 0.35 0.35]; %need to base these on real data

%%

useConstantWVel = 0;

%set melt parameters
casename = {'Seasonal'};

%relative velocity
load('./Infiles/ADCP/ADCP_seasonal.mat', 'zadcp','vadcp','tadcp','wvel'); %use Becca's summer exchange flow

%wave, turb water, turb air, free air, free water
calcMelt = [1 1 1 1 1];

%% 

dz = 5;

L = 50:50:1000; %iceberg size classes (50 m increments)

ni = length(L);

%% 

startDate = datenum(2008,1,1,0,0,0);

modelDays = 365;
timeSpan = 86400 .* modelDays; %number of days to run simulation

interpTime = startDate:1:startDate + modelDays;

[monthlyTime it] = sort(tdmon);

%%

oceanTemp = Tseas_mon(:,it); %ocean temperature (deg C)
oceanSal = Sseas_mon(:,it); %ocean salinity
oceanTSDepth = ctdz; %ocean depth (m)

airTemp = MONmat(it); %air temp (deg C)
swFlux = MONmsw(it); %shortwave radiation flux (W/m^-2
windSpeed = MONmspda(it); %wind speed (m s^-1)
IceConc = MONcice(it); %sea-ice concentration (0-1)

oceanWVel = wvel;

if useConstantWVel

    oceanWVel = (oceanWVel .* 0) + 0.03;

end

oceanVelDepth = zadcp;

oceanVelTime = datenum(2008,1,1,0,0,0):1:(datenum(2009,1,1,0,0,0)-1);

[year month day hours minutes seconds] = datevec(oceanVelTime);

uniqueMonth = unique(month);

oceanUVel = [];

for i = 1:length(uniqueMonth)
    
    im = find(month == uniqueMonth(i));
    
    field = vadcp(:,1:length(im),i);
    
    oceanUVel = cat(2,oceanUVel,field);
    
end

%% 
%interp onto daily timestep

forcing.time = interpTime;

forcing.airTemp = interp1(monthlyTime,airTemp,interpTime);
forcing.swFlux = interp1(monthlyTime,swFlux,interpTime);
forcing.windSpeed = interp1(monthlyTime,windSpeed,interpTime);
forcing.iceConc = interp1(monthlyTime,IceConc,interpTime);

forcing.oceanTemp = interp2(monthlyTime,oceanTSDepth,oceanTemp,interpTime,oceanVelDepth);
forcing.oceanSal = interp2(monthlyTime,oceanTSDepth,oceanSal,interpTime,oceanVelDepth);

forcing.oceanUVel = interp2(oceanVelTime,oceanVelDepth,oceanUVel,interpTime,oceanVelDepth);

forcing.oceanWVel = interp1(monthlyTime,oceanWVel,interpTime);
forcing.oceanWVel = repmat(forcing.oceanWVel,[length(oceanVelDepth) 1]);

forcing.oceanTSDepth = oceanVelDepth;
forcing.oceanVelDepth = oceanVelDepth;

forcing.speckleFe = speckleFe;

useRealFeConc = 1;

%% 

%Sutherland et al. 2014 residence times
%81 ± 67 
%meanRes = 81;
%stdRes = 67;

plotFigure = 0;
saveMovie = 0;

startTime = datenum(2008,5,1,0,0,0);
is = find(forcing.time >= startTime,1);
ie = 349;

runDays = length(is:ie);
forcing.is = is;

forcing.airTemp = forcing.airTemp(is:ie);
forcing.swFlux = forcing.swFlux(is:ie);
forcing.windSpeed = forcing.windSpeed(is:ie);
forcing.iceConc = forcing.iceConc(is:ie);
forcing.oceanTemp = forcing.oceanTemp(:,is:ie);
forcing.oceanSal = forcing.oceanSal(:,is:ie);
forcing.oceanUVel = forcing.oceanUVel(:,is:ie);
forcing.oceanWVel = forcing.oceanWVel(:,is:ie);

for i = 2:2
   
    icebergType = i; %1 = shell, 2 = speckle, 3 = basal layer only
    
    for j = 1:length(L)
        
        timeSpan = runDays.* 86400;
        
        disp(['Length: ' num2str(L(j))]);
        
        iceberg(j) = iceberg_melt(L(j),dz,timeSpan,forcing,calcMelt,icebergType,useRealFeConc,plotFigure,saveMovie);
        
    end
    
    save([saveDir 'iceberg_' num2str(i) 'wave_melt.mat'],'-v7.3');
    
    clear iceberg
    
end

%% 
