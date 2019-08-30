function out = iceberg_melt(L,dz,timeSpan,forcing,calcMelt,icebergType,useRealFeConc,plotFigure,saveMovie)

matrixSaveDir = '/Users/carrolld/Documents/research/greenland/iceberg_iron/mat/moon_iceberg_model/iceberg_matrix/';
movieDir = '/Users/carrolld/Documents/research/greenland/iceberg_iron/movies/';

%homogenousFe = 9.274 .* 10^-6; %9.274 umol Fe

if useRealFeConc
    
    shellFeConc = 100 .* 10^-6; %100 umol Fe
    interiorFeConc = 167 .* 10^-9; %167 nmol Fe
    
else
    
    shellFeConc = 2;
    interiorFeConc = 1;
    
end

netFeFrac = 0.09;

useShell = 0;
useSpeckle = 0;
useBasal = 0;

if icebergType == 1
    
    useShell = 1;
    
elseif icebergType == 2
    
    useSpeckle = 1;
    
elseif icebergType == 3
    
    useBasal = 1;
    
end

useIntegerIceShape = 1; %round initial iceberg dimensions

calcWaveMelt  = 0; %calcMelt(1); %wave erosion
calcWaveMelt  = calcMelt(1); %wave erosion

calcTurbWaterMelt = calcMelt(2); %forced underwater convection, acts on iceberg side and base
calcTurbAirMelt = calcMelt(3); %forced convection in air, acts on iceberg sides and top
calcFreeAirMelt = calcMelt(4); %melting in air, reduces thickness only
calcFreeWaterMelt = calcMelt(5); %buoyant convection, only on iceberg sides

%%

padX = 30;

ice_init = init_iceberg_size(L,dz);  %initialize iceberg sizes and geometry

%round iceberg dimensions to avoid fractions in initial iceberg matrix
if useIntegerIceShape
    
    ice_init.L = round2even(ice_init.L);
    ice_init.W = round2even(ice_init.W);
    
    ice_init.freeB = round2even(ice_init.freeB);
    ice_init.K = round2even(ice_init.K);
    ice_init.TH = ice_init.freeB + ice_init.K;
    
    for i = 1:length(ice_init.uwL)
        
        ice_init.uwL(i) = round2even(ice_init.uwL(i));
        
    end
    
    for i = 1:length(ice_init.uwW)
        
        ice_init.uwW(i) = round2even(ice_init.uwW(i));
        
    end
    
    maxL = round2even(nanmax([ice_init.L; ice_init.uwL])) + padX;
    maxW = round2even(nanmax([ice_init.W; ice_init.uwW])) + padX;
    
else
    
    %only round max L/W matrix dimensions to even number
    maxL = round2even(nanmax([ice_init.L; ice_init.uwL])) + padX;
    maxW = round2even(nanmax([ice_init.W; ice_init.uwW])) + padX;
    
end

maxFreeB = ice_init.freeB+padX;
maxThickness = ice_init.TH + padX;

ice_init.maxL = maxL;
ice_init.maxW = maxW;
ice_init.maxThickness = maxThickness;
ice_init.maxFreeB = maxFreeB;
ice_init.even = useIntegerIceShape;

[iceRunFrac iceInitVol iceInitSA interpInitProfile] = generate_iceberg_matrix(ice_init,[],[]);

feVol = iceInitVol .* netFeFrac;

%%

%dz = 1;

nz = 500;

dt = 86400; %all melt rates in m s^-1

t = dt:dt:timeSpan;

t = t'; %time vector from 0 to total timeSpan in seconds

nt = length(t); %number of timesteps, first timestep is 0

%%

ni = 1;

Mwave  = zeros(ni,nt); %melt volume for waves, affects just top layer

mw = Mwave;
ma = mw;
ms = mw;

Mturbw = zeros(nz,nt); %forced convection underwater, acts on side and base
Mturba = zeros(nt); %forced convection in air, acts on sides and top
Mfreea = zeros(nt); %melting in air, reduces thickness only
Mfreew = zeros(nz,nt); %buoyant convection, only on sides

mtw = Mturbw;
mb = mtw;

%%
%ocean velocity, horizontal currents and vertical velocities (plumes)

%find keel depth
ik = dsearchn(forcing.oceanVelDepth,ceil(ice_init.K));

%mean velocity
forcing.meanUVel = repmat(nanmean(forcing.oceanUVel),[length(forcing.oceanVelDepth) 1]);

%remove mean
forcing.oceanUVel = forcing.oceanUVel - forcing.meanUVel;

%add in vertical velocity
forcing.oceanUVel = forcing.oceanUVel + forcing.oceanWVel;

%speed
forcing.oceanUVel = abs(forcing.oceanUVel);

%%

shellThickness = (iceInitVol .* netFeFrac) ./ iceInitSA; %thickness of debris shell (m)
%shellThickness = 5;

scale = 2;

if useShell
    
    ice_fe_init.L = ice_init.L - (shellThickness .* scale);
    ice_fe_init.W = ice_init.W - (shellThickness .* scale);
    ice_fe_init.uwL = ice_init.uwL - (shellThickness .* scale);
    ice_fe_init.uwW = ice_init.uwW - (shellThickness .* scale);
    
    ice_fe_init.interp_uwL = interpInitProfile.uwL - (shellThickness .* scale);
    ice_fe_init.interp_uwW = interpInitProfile.uwW - (shellThickness .* scale);
    
    ice_fe_init.interp_uwL(ice_fe_init.interp_uwL < 0) = 0;
    ice_fe_init.interp_uwW(ice_fe_init.interp_uwW < 0) = 0;
    
    ice_fe_init.Z = ice_init.Z;
    ice_fe_init.freeB = ice_init.freeB - (shellThickness);
    ice_fe_init.K = ice_init.K - (shellThickness);
    ice_fe_init.TH = (ice_fe_init.freeB + ice_fe_init.K);
    
end

if useSpeckle
    
    ice_fe_init.L = ice_init.L;
    ice_fe_init.W = ice_init.W;
    ice_fe_init.uwL = ice_init.uwL;
    ice_fe_init.uwW = ice_init.uwW;
    
    ice_fe_init.interp_uwL = interpInitProfile.uwL;
    ice_fe_init.interp_uwW = interpInitProfile.uwW;
    
    ice_fe_init.interp_uwL(ice_fe_init.interp_uwL < 0) = 0;
    ice_fe_init.interp_uwW(ice_fe_init.interp_uwW < 0) = 0;
    
    ice_fe_init.Z = ice_init.Z;
    ice_fe_init.freeB = ice_init.freeB;
    ice_fe_init.K = ice_init.K;
    ice_fe_init.TH = (ice_fe_init.freeB + ice_fe_init.K);
    
end

if useBasal
    
    ice_fe_init.L = ice_init.L;
    ice_fe_init.W = ice_init.W;
    ice_fe_init.uwL = ice_init.uwL;
    ice_fe_init.uwW = ice_init.uwW;
    ice_fe_init.Z = ice_init.Z;
    ice_fe_init.freeB = ice_init.freeB;
    
    ice_fe_init.interp_uwL = interpInitProfile.uwL;
    ice_fe_init.interp_uwW = interpInitProfile.uwW;
    
    ice_fe_init.interp_uwL(ice_fe_init.interp_uwL < 0) = 0;
    ice_fe_init.interp_uwW(ice_fe_init.interp_uwW < 0) = 0;
    
    interpArea = (ice_fe_init.interp_uwL .* ice_fe_init.interp_uwW);
    
    interpArea = cumsum(fliplr(interpArea));
    
    maxUwVol = max(interpArea);
    %maxTotalVol = maxUwVol + (ice_fe_init.L .* ice_fe_init.W .* ice_fe_init.freeB);
    maxTotalVol = feVol;
    
    interpArea(interpArea == 0) = [];
    
    iz = find(interpArea > (netFeFrac.*maxTotalVol),1); %depth index
    
    iiz = (length(ice_fe_init.interp_uwL) - iz); %index from basal layer
    
    basalInt = interpArea(iz);
    basalFrac = feVol - basalInt;
    
    basalThicknessInt = iz;
    basalThicknessFrac = basalFrac ./ (interpArea(iz)+1);
    basalThickness = basalThicknessInt + basalThicknessFrac;
    
    ice_fe_init.K = ice_init.K - basalThickness;
    
    ice_fe_init.TH = (ice_fe_init.freeB + ice_fe_init.K);
    
end

ice_fe_init.maxL = maxL;
ice_fe_init.maxW = maxW;
ice_fe_init.maxThickness = maxThickness;
ice_fe_init.maxFreeB = maxFreeB;

ice_fe_init.even = 0; %useIntegerIceShape;

[iceFeConc iceFeVol iceFeSA interpFeProfile] = generate_iceberg_matrix(ice_fe_init,ice_fe_init.interp_uwL,ice_fe_init.interp_uwW);

%iceFeConc(:,:,padX) = 1;

%%

clear field

iceShellDiff = iceRunFrac - iceFeConc;

shellIndex = find(iceShellDiff ~= 0);

if useShell
    
    field = iceRunFrac;
    
    wetMask = find(iceFeConc == 0);
    field(wetMask) = nan;
    
    interiorIndex = find(~isnan(field));
    field(interiorIndex) = interiorFeConc;
    
    field(shellIndex) = shellFeConc .* iceShellDiff(shellIndex);
    
end

if useSpeckle
    
    field = iceRunFrac;
    
    wetMask = find(iceRunFrac == 0);
    field(wetMask) = nan;
    
    speckleIndex1 = find(~isnan(field));
    field(speckleIndex1) = interiorFeConc;
    
    speckleIndex2 = datasample(speckleIndex1,10^9);
    
    randIceRunFrac = iceRunFrac(speckleIndex2);
    randSpeckleFeConc = datasample(forcing.speckleFe,length(speckleIndex2));
    
    tempFeConc = cumsum(randSpeckleFeConc .* randIceRunFrac);
    
    ix = find(tempFeConc >=  (shellFeConc .* feVol),1);
    
    field(speckleIndex2(1:ix)) = randSpeckleFeConc(1:ix);
    
    clear wetMask speckleIndex 1 speckleIndex2 randIceRunFrac randSpeckleFeConc tempFeConc ix 
    
end

if useBasal
    
    field = iceRunFrac;
    
    wetMask = find(iceFeConc == 0);
    field(wetMask) = nan;
    
    interiorIndex = find(~isnan(field));
    field(interiorIndex) = interiorFeConc;
    
    field(shellIndex) = shellFeConc .* iceShellDiff(shellIndex);
    field(:,:,iiz+1) = field(:,:,iiz+1) + interiorFeConc;
    
end

%%

iceFeConc = field;

%set bounds of iceberg matrix to initial shape
ice_run.maxL = maxL;
ice_run.maxW = maxW;
ice_run.maxThickness = maxThickness;
ice_run.maxFreeB = maxFreeB;
ice_run.even = 0; %useIntegerIceShape;

clear field

%%

folderName = ['iceberg_L' num2str(L) '_type_' num2str(icebergType)];

ice = ice_init;
use ice %brings structure variables into top level of workspace

dz = 1;
Z = 1:500;

uwL = interpInitProfile.uwL;
uwW = interpInitProfile.uwW;

iceRunFeVol(1) = nansum(nansum(nansum(iceRunFrac .* iceFeConc)));

[m n] = size(squeeze(iceRunFrac(maxW/2,:,:)));
[xx zz] = meshgrid(1:m,1:n);
zz = -zz + ice_init.freeB + padX + 1;

if plotFigure

    mkdir([movieDir folderName]);

    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    
end

for j = 1:nt %time loop
    
    %start calculating melt, get melt rates for each process including, then update at end
    
    iKeel = ceil(ice.K ./ dz); % index of deepest iceberg level
    
    if (calcWaveMelt)
        
        SST = nanmean(forcing.oceanTemp(1:10,j,1)); %0-10 m
        
        WH(j) = 0.010125 .* ((abs(forcing.windSpeed(j))) .^ 2); %assume wind vel >> ocean vel, this estimates wave height
        WHdep = min(freeB, 1 .* WH(j)); %thickness that wave height affects
        
        %apply 1/2 mw to L and 1/2 mw to uwL(1,:)
        mw(j) = melt_wave(forcing.windSpeed(j),SST,forcing.iceConc(j)); %already in m s^-1
        
        mw(j) = mw(j) .* dt; %melt rate in m s^-1, this now in m day^-1
        lint = nanmean([L uwL(1)]); %mean of length and first layer underwater
        
        Mwave(j) = 1.*(mw(j) .* WHdep .* lint)   ... %1 lengths
            + 1.*(mw(j) .* WHdep .* lint); %+ 1 widths (coming at it obliquely)
        
        %base on wave height estimate, to get right volume taken off-- but doesn't do L right then! FIX
        mwabove = WHdep ./ freeB; %ratio of wave height to freeB
        mwbelow = WHdep ./ dz; %ratio of wave height to first layer thickness
        
    else
        
        mw(j) = 0;
        mwabove = 0;
        mwbelow = 0;
        
    end
    
    if (calcTurbWaterMelt)
        
        %do for each depth level where iceberg exists
        for k = 1:iKeel
            
            T_far = interp1(forcing.oceanTSDepth,forcing.oceanTemp(:,j),Z(k));
            S_far = interp1(forcing.oceanTSDepth,forcing.oceanSal(:,j),Z(k));
            U_rel = interp1(forcing.oceanVelDepth,forcing.oceanUVel(:,j),Z(k));
            
            [mtw(k,j), ~, ~] = melt_forcedwater(T_far,S_far,Z(k),U_rel);
            mtw(k,j) = mtw(k,j) .* dt;
            
            Mturbw(k,j) = 2.*(mtw(k,j) .* dz .* uwL(k)) ... %2 lengths
                + 1.*(mtw(k,j) .* dz .* uwW(k)); %+ 1 widths (don't count lee side)
            
        end
        
        U_rel_keel = interp1(forcing.oceanVelDepth,forcing.oceanUVel(:,j),Z(iKeel));
        
        [mtw(iKeel,j), ~, ~] = melt_forcedwater(T_far,S_far,Z(iKeel),U_rel_keel);
        mtw(iKeel,j) = mtw(iKeel,j) .* dt;
        
        dz_keel = -((iKeel-1).*dz - K); %final layer depth
        
        %do keel layer now
        Mturbw(iKeel,j) = 2.*(mtw(iKeel,j) .* dz_keel .* uwL(iKeel)) ... %2 lengths
            + 1.*(mtw(iKeel,j) .* dz_keel .* uwW(iKeel)); %+ 1 widths (don't count lee side)
        
    else
        
        mtw(1:nz,j) = 0;
        
    end
    
    if (calcTurbAirMelt)
        
        ma(j) = melt_forcedair(forcing.airTemp(j),forcing.windSpeed(j),L);
        ma(j) = ma(j) .* dt; %melt rate in m s^-1
        
        Mturba(j) = 2 .* (ma(j) .* dz .* L) ... %2 lengths
            + 1 .* (ma(j) .* dz .* W) ... %+ 1 widths (don't count lee side)
            + 0.5 .* (ma(j) .* L .* W); %+ 1/2 top surface area (account for not flat)
        
    else
        
        ma(j) = 0;
        
    end
    
    if (calcFreeAirMelt)
        
        %determine debris level on top layer of freeboard
        
        surfFlag = 0;
        
        c = 1;
        
        while(surfFlag ~= 1)
            
            top = iceFeConc(:,:,c);
            
            top(isnan(top)) = [];
            
            temp = unique(top);
            
            %if above freeboard in padded grid region
            if(isempty(temp))
                
                temp  = 0;
                
            end
            
            %if empty layer is found, iterate downwards
            if(temp == 0)
                
                surfFlag = 0;
                clear temp
                
                c = c + 1;
                
            %if solid shell
            elseif(temp == shellFeConc)
                
                surfFlag = 1;
                %alb = 0.3; %albedo of dirty ice
                alb = 0.7; %albedo of clean ice
                
                %if speckled / interior
            else
                
                surfFlag = 1;
                alb = 0.7; %albedo of clean ice
                
            end
            
        end
        
        clear temp
        
        ms(j) = melt_solar(forcing.swFlux(j),alb);
        
        ms(j) = ms(j) .* dt; %melt rate in m s^-1
        
        Mfreea(j) = (ms(j) .* W .* L); %acts on top surface area only
        
    else
        
        ms(j) = 0;
        
    end
    
    if (calcFreeWaterMelt)
        
        %do for each depth level where iceberg exists
        for k = 1:iKeel-1
            
            T_far = interp1(forcing.oceanTSDepth,forcing.oceanTemp(:,j),Z(k));
            S_far = interp1(forcing.oceanTSDepth,forcing.oceanSal(:,j),Z(k));
            
            mb(k,j) = melt_buoyantwater(T_far,S_far,'cis'); %using Bigg method, then S doesn't matter
            mb(k,j) = mb(k,j) .* dt;
            
            Mfreew(k,j) = 2 .* (mb(k,j) .* dz .* uwL(k)) ... %2 lengths
                + 2 .* (mb(k,j) .* dz .* uwW(k)); %+ 2 widths
            
        end
        
        %do keel layer
        dz_keel = -((iKeel-1).*dz - K); %final layer depth
        
        Mfreew(iKeel,j) = 2 .* (mb(iKeel,j) .* dz_keel .* uwL(iKeel)) ... %2 lengths
            + 2 .* (mb(iKeel,j) .* dz_keel .* uwW(iKeel)); %+ 2 widths
        
    else
        
        mb(1:nz,j) = 0;
        
    end
    
    %reduce thickness from top and bottom
    freeB = freeB - ms(j) - ma(j);
    K = K - mtw(iKeel,j);
    TH = K + freeB;
    
    %reduce thickness on sides, do one L and update W's accordingly
     
    %uwL(1) = uwL(1) - (scale.* mw(j) .* (mwbelow / 1)); %length under waterline
    %uwW(1) = uwW(1) - (scale.* mw(j) .* (mwbelow / 1)); %length under waterline
    
    %putting all mw at L means taking out too much volume, b/c it is freeB high
    %L = L - (scale .* ma(j)) - (scale .* mw(j) .* (mwabove / 1)); %length at waterline
    
    %turbulent and free convection melt at depth
    %for k=1:iKeel
        
    %    uwL(k) = uwL(k) - (scale .* mtw(k,j)) - (scale .* mb(k,j));
    %    uwW(k) = uwW(k) - (scale .* mtw(k,j)) - (scale .* mb(k,j));
        
    %end

    L = L - (scale .* ma(j)) - (scale .* mw(j));
    
    %turbulent and free convection melt at depth
    for k=1:iKeel
        
        uwL(k) = uwL(k) - (scale .* mtw(k,j)) - (scale .* mb(k,j)) - (scale .* mw(j));
        uwW(k) = uwW(k) - (scale .* mtw(k,j)) - (scale .* mb(k,j)) - (scale .* mw(j));
        
    end

    %%
    %****FIX!****
    
    uwW = uwL ./ 1.62; %update widths accordingly
    W = L ./ 1.62; %update widths
    
    %%
    
    %update geometry
    rhoi = 900;
    rati = rhoi / 1024; %ratio of ice density to water density
    
    kindnew = ceil(K ./ dz); %check to see if it got rid of level
    
    if (kindnew < iKeel)  % then get rid of that bottom layer
        
        uwL(iKeel) = nan;
        uwW(iKeel) = nan;
        uwV(iKeel) = nan;
        iKeel = kindnew;
        
    end
    
    uwV(1:(iKeel-1)) = dz .* uwL(1:(iKeel-1)) .* uwW(1:(iKeel-1));
    dzk(j) = -((iKeel-1).*dz - K);
    uwV(iKeel) = dzk(j) .* uwL(iKeel) .* uwW(iKeel); %bottom layer is not dz thick
    
    sailV = freeB .* L .* W;
    totalV = nansum(uwV) + sailV;
    
    sailV = (1 - rati) .* totalV;
    
    freeB = sailV ./ (L .* W); %new freeboard

    K = TH - freeB;
    
    wstab = 0.92;
    lwrat = L ./ TH;  % length to thickness ratio
    
    %disp(num2str(lwrat));
    
    if (lwrat < wstab || uwL(1) <= 0)  % if L/thickness < 0.7, then roll
        
        %disp('unstable');
        
        out.ice_init = ice_init;
        
        %integrated melt terms
        out.Mwave = Mwave;
        out.Mfreea = Mfreea;
        out.Mturbw = Mturbw;
        out.Mturba = Mturba;
        out.Mfreew = Mfreew;
        
        %melt terms
        out.i_mwave = mw; %m day^-1
        out.i_mfreea = ms; %m day^-1
        out.i_mturbw = mtw; %m day^-1
        out.i_mturba = ma; %m day^-1
        out.i_mfreew = mb; %m day^-1
        
        %mean over all time, depths, processes in m day^-1
        out.i_mtotalm = nanmean(mw(:)) + nanmean(mb(:)) + nanmean(ms(:)) + ...
            nanmean(ma(:)) + nanmean(mtw(:));
        
        out.iceRunVol = iceRunVol;
        out.iceRunFeVol = iceRunFeVol;
        
        out.t = t;
        
        out.meltparams.calcWaveMelt  = calcWaveMelt;
        out.meltparams.calcTurbWaterMelt  = calcTurbWaterMelt;
        out.meltparams.calcTurbAirMelt  = calcTurbAirMelt;
        out.meltparams.calcFreeAirMelt= calcFreeAirMelt;
        out.meltparams.calcFreeWaterMelt = calcFreeWaterMelt;
        
        out.meltparams.Urel = forcing.oceanUVel;
        out.meltparams.ctdtemp = forcing.oceanTemp;
        out.meltparams.ctdsalt = forcing.oceanSal;
        out.meltparams.ctdz = forcing.oceanTSDepth;
        out.meltparams.windvel = forcing.windSpeed;
        out.meltparams.seaiceConc = forcing.iceConc;
        out.meltparams.airTemp = forcing.airTemp;
        out.meltparams.solarins = forcing.swFlux;
        
        return
        
    end
    
    %update iceberg matrix geometry
    ice_run.L = L;
    ice_run.W = W;
    ice_run.TH = TH;
    ice_run.uwL = uwL;
    ice_run.uwW = uwW;
    ice_run.K = K;
    ice_run.Z = Z;
    ice_run.freeB = freeB;
    ice_run.even = 0; %useIntegerIceShape;
    
    uwL(uwL < 0) = 0;
    uwW(uwW < 0) = 0;
    
    [iceRunFrac IceRunVol iceSA interpRunProfile] = generate_iceberg_matrix(ice_run,uwL,uwW);
        
    iceRunFe = iceRunFrac .* iceFeConc;

    save([matrixSaveDir 'iceberg_matrix_' num2str(j)],'iceRunFe');
    
    iceRunVol(j) = nansum(nansum(nansum(iceRunFrac)));
    iceRunFeVol(j) = nansum(nansum(nansum(iceRunFe)));
    
    %%
    
    if plotFigure
        
        fs = 30;
        lw = 3;
        bgColor = [85 174 255] / 255;
        
        cc1 = subplot(121);
        
        hold on
        
        set(gca,'color',bgColor);
        
        temp = squeeze(iceRunFrac(maxL/2,:,:));
        temp(temp == 0) = nan;
        pcolor(xx,zz,temp');
        shading flat
        
        plot(xx(1,:),(xx(1,:) .* 0),'LineWidth',lw','Color','k','LineStyle','--');
        
        colormap(cc1,flipud(cbrewer('seq','Blues',500)));
        
        hcb1 = colorbar;
        set(get(hcb1,'ylabel'),'String','Ice Cell Volume (m^-^3)','FontSIze',fs);
        
        caxis([0 1]);
        
        xlim([min(min(xx)) max(max(xx))+1]);
        ylim([min(min(zz)) max(max(zz))]);
        
        xlabel('Distance (m)');
        ylabel('Height (m)');
        
        box on
        set(gca,'GridLineStyle','--');
        set(gca,'LineWidth',2);
        set(gca,'FontSize',fs);
        
        title(['Initial Length = ',num2str(ice_init.L) ', Day: ' num2str(j)],'FontWeight','Normal');
        
        cc2 = subplot(122);
        
        hold on
        
        set(gca,'color',bgColor);
        
        temp = squeeze(iceFeConc(maxL/2,:,:) .* iceRunFrac(maxL/2,:,:));
        temp2 = temp;
        it = find(isnan(temp2));
        temp2(it) = 0;
        temp2 = real(log10(temp2));
        temp2(isinf(temp2)) = nan;
        
        if useRealFeConc
            
            pcolor(xx,zz,temp2');
            
        else
            
            temp(temp == 0) = nan;
            pcolor(xx,zz,temp');
            
        end
        
        shading flat
        
        plot(xx(1,:),(xx(1,:) .* 0),'LineWidth',lw','Color','k','LineStyle','--');
        
        colormap(cc2,flipud(cmocean('grey',500)));
        
        hcb2 = colorbar;
        set(get(hcb2,'ylabel'),'String','log10(Fe) (mol Fe m^-^3)','FontSize',fs);
        
        if useRealFeConc
            
            caxis([-9 -4]);
            
        else
            
            caxis([0 2]);
            
        end
        
        xlim([min(min(xx)) max(max(xx))+1]);
        ylim([min(min(zz)) max(max(zz))]);
        
        xlabel('Distance (m)');
        
        box on
        set(gca,'GridLineStyle','--');
        set(gca,'LineWidth',2);
        set(gca,'FontSize',fs);
        
        drawnow
        
        if saveMovie
            
            export_fig(hFig1,[movieDir folderName '/figure_ ' num2str(j) '.png'],'-png'); %, '-r300');
            
        end
        
        clear temp temp2
        
        %clf
        
    else
        
        hold on
        %plot(uwL,'k')
        plot(interpRunProfile.uwL,'r')
        
        drawnow
        
    end
    
    %pause
    clf
    
    clear iceRunFrac iceRunFe
    %disp(num2str(j));
    
end

%%
%output

out.ice_init = ice_init;

%integrated melt terms
out.Mwave = Mwave;
out.Mfreea = Mfreea;
out.Mturbw = Mturbw;
out.Mturba = Mturba;
out.Mfreew = Mfreew;

%melt terms
out.i_mwave = mw; %m day^-1
out.i_mfreea = ms; %m day^-1
out.i_mturbw = mtw; %m day^-1
out.i_mturba = ma; %m day^-1
out.i_mfreew = mb; %m day^-1

%mean over all time, depths, processes in m day^-1
out.i_mtotalm = nanmean(mw(:)) + nanmean(mb(:)) + nanmean(ms(:)) + ...
    nanmean(ma(:)) + nanmean(mtw(:));

out.iceRunVol = iceRunVol;
out.iceRunFeVol = iceRunFeVol;

out.t = t;

out.meltparams.calcWaveMelt  = calcWaveMelt;
out.meltparams.calcTurbWaterMelt  = calcTurbWaterMelt;
out.meltparams.calcTurbAirMelt  = calcTurbAirMelt;
out.meltparams.calcFreeAirMelt= calcFreeAirMelt;
out.meltparams.calcFreeWaterMelt = calcFreeWaterMelt;

out.meltparams.Urel = forcing.oceanUVel;
out.meltparams.ctdtemp = forcing.oceanTemp;
out.meltparams.ctdsalt = forcing.oceanSal;
out.meltparams.ctdz = forcing.oceanTSDepth;
out.meltparams.windvel = forcing.windSpeed;
out.meltparams.seaiceConc = forcing.iceConc;
out.meltparams.airTemp = forcing.airTemp;
out.meltparams.solarins = forcing.swFlux;

%%
%clean up

clear i j k ni nz nt kind wstab lwrat rati do_*
