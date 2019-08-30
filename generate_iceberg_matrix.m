function [iceFrac iceVol iceSA interpProfile] = generate_iceberg_matrix(ice_init,uwL,uwW)

%compute freeboard matrix

freeLInt = floor(ice_init.L);
freeLFrac = freeLInt - ice_init.L;

freeWInt = floor(ice_init.W);
freeWFrac = freeWInt - ice_init.W;

freeBInt = floor(ice_init.freeB);
freeBFrac = abs(freeBInt - ice_init.freeB);

freeLIntSide = freeLInt / 2;
freeLIntFracSide = (floor(freeLIntSide) - freeLIntSide) + freeLFrac/2;
freeLIntSide = floor(freeLIntSide);

freeWIntSide = freeWInt / 2;
freeWIntFracSide = (floor(freeWIntSide) - freeWIntSide) + freeWFrac/2;
freeWIntSide = floor(freeWIntSide);

freeX = [[abs(freeLIntFracSide) ones(1,freeLIntSide)]  fliplr([abs(freeLIntFracSide) ones(1,freeLIntSide)])];
freeY = [[abs(freeWIntFracSide) ones(1,freeWIntSide)] fliplr([abs(freeWIntFracSide) ones(1,freeWIntSide)])];

if freeBFrac == 0
    
    freeZ = [1 ones(1,freeBInt)];

else
        
    freeZ = [freeBFrac ones(1,freeBInt)];

end

padX = ceil((ice_init.maxL - length(freeX)) / 2);
padY = ceil((ice_init.maxW - length(freeY)) / 2);
padZ = ceil((ice_init.maxFreeB - length(freeZ)));

freeZ = [zeros(1,padZ) freeZ];

[xx yy] = meshgrid([zeros(1,padX) freeX zeros(1,padX)],[zeros(1,padY) freeY zeros(1,padY)]);

%make 2D matrix of ice conc
free2D = xx .* yy;

clear xx yy

[m n] = size(free2D);

free3D = repmat(free2D,[1 1 length(freeZ)]);

freeZZ = repmat(freeZ',[1 m n]);
freeZZ = permute(freeZZ,[2 3 1]);

freeIceFrac = free3D .* freeZZ;

%% 

dz = 5;

%subsurface

subThickInt = floor(ice_init.K);
subThickFrac = abs(subThickInt - ice_init.K);

%number of z-cells to interpolate
numSubZ = ceil(subThickFrac) + subThickInt;

if(~isempty(uwL) && ~isempty(uwW))
    
    interpSub_uwL = uwL;
    interpSub_uwW = uwW;
    
    interpSub_uwL(subThickInt+1:end) = 0; 
    interpSub_uwW(subThickInt+1:end) = 0; 
    
    interpSub_uwL(interpSub_uwL < 0) = 0;
    interpSub_uwW(interpSub_uwW < 0) = 0;
    
else
    
    startDepth = 1;
    
    %interpolate down to initial iceberg thickness
    interpSubZ = startDepth:ice_init.maxThickness;
    
    interpSub_uwL =  interp1(ice_init.Z,ice_init.uwL,interpSubZ);
    interpSub_uwW = interp1(ice_init.Z,ice_init.uwW,interpSubZ);
    
    interpSub_uwL(isnan(interpSub_uwL)) = 0;
    interpSub_uwW(isnan(interpSub_uwW)) = 0;
    
    %fill < dz
    interpSub_uwL(1:dz-1) = interpSub_uwL(dz);
    interpSub_uwW(1:dz-1) = interpSub_uwW(dz);
    
    interpSub_uwL(numSubZ+1:end) = 0;
    interpSub_uwW(numSubZ+1:end) = 0;
    
end

ik = find(interpSub_uwL == 0,1);

if (ice_init.even == 1)
    
    for i = 1:length(interpSub_uwL)
        
        interpSub_uwL(i) = round2even(interpSub_uwL(i));
        
    end
    
    for i = 1:length(interpSub_uwW)
        
        interpSub_uwW(i) = round2even(interpSub_uwW(i));
        
    end
    
end

for i = 1:ice_init.maxThickness
    
    subLInt = floor(interpSub_uwL(i));
    subLFrac = subLInt - interpSub_uwL(i);
 
    subWInt = floor(interpSub_uwW(i));
    subWFrac = subWInt - interpSub_uwW(i);
    
    subLIntSide = subLInt ./ 2;
    subLIntFracSide = (floor(subLIntSide) - subLIntSide) + (subLFrac ./ 2);
    subLIntSide = floor(subLIntSide);

    subWIntSide = subWInt ./ 2;
    subWIntFracSide = (floor(subWIntSide) - subWIntSide) + (subWFrac ./ 2);
    subWIntSide = floor(subWIntSide);

    subX = [[abs(subLIntFracSide) ones(1,subLIntSide)]  fliplr([abs(subLIntFracSide) ones(1,subLIntSide)])];
    subY = [[abs(subWIntFracSide) ones(1,subWIntSide)] fliplr([abs(subWIntFracSide) ones(1,subWIntSide)])];
    
    padX = (n - length(subX)) / 2;
    padY = (m - length(subY)) / 2;
    
    [xx yy] = meshgrid([zeros(1,padX) subX zeros(1,padX)],[zeros(1,padY) subY zeros(1,padY)]);
    
    sub2D = xx .* yy;
    
    subIceFrac(:,:,i) = sub2D;
    
    clear subX subY padX padY sub2D
    
end

if ik > 1

    subIceFrac(:,:,ik-1) = subIceFrac(:,:,ik-1) .* subThickFrac;
    iceSubBaseSA = interpSub_uwL(ik-1) .* interpSub_uwW(ik-1);

else
    
    iceSubBaseSA = 0;
    
end
    
iceFrac = cat(3,freeIceFrac,subIceFrac);
%iceFrac(iceFrac == 0) = nan;

iceVol =  nansum(nansum(nansum(freeIceFrac))) + nansum(nansum(nansum(subIceFrac)));

%compute surface area
iceFreeSA = (((ice_init.L .* 2) + (ice_init.W .* 2)) .* freeBInt) + (ice_init.L .* ice_init.W);
iceSubSA = nansum(((interpSub_uwL .* 2) + (interpSub_uwW .* 2)));

iceSA = iceFreeSA + iceSubSA + iceSubBaseSA;

interpProfile.uwL = interpSub_uwL;
interpProfile.uwW = interpSub_uwW;

clear freeIceFrac suIceFrac
    
end
