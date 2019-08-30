function icebergs = barker_carea(L, K, dz, LWratio)
%
% calculates underwater cross sectional areas using Barker et al. 2004,
% and converts to length underwater (for 10 m thicknesses, this is just CrossArea/10) and sail area for K<200, 
% for icebergs K>200, assumes tabular shape
%
% [CArea, UWlength, SailA] = barker_carea(L)
%
% L is vector of iceberg lengths, 
% K is keel depth
% (if narg<2), then it calculates keel depths from this using keeldepth.m
% dz = layer thickness to use
% LWratio: optional argument, default is 1.62:1 L:W, or specify
%
% all variables in structure "icebergs"
%   CA is cross sectional area of each 10 m layer underwater
%   uwL is length of this 10 m underwater layer
%   Z is depth of layer
%   uwW calculated from length using length to width ratio of 1.62:1
% 
%   also get volumes and masses
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first get keel depth
L = L(:);

if nargin == 1

    K = keeldepth(L,'barker'); % K = keeldepth(L,'mean');
    dz = 10; 
    LWratio = 1.62;

elseif nargin<4

    LWratio = 1.62; 

end

K = K(:);

% barker table 5 (0-200 m)
if(dz == 10) % originally for dz=10 m layers
    
    a = [9.51;11.17;12.48;13.6;14.3;13.7;13.5;15.8;14.7;11.8;11.4;10.9;10.5;10.1;9.7;9.3;8.96;8.6;8.3;7.95];
    b = -[25.9;107.5;232;344.6;457;433;520;1112;1125;853;931;1007;1080;1149;1216;1281;1343;1403;1460;1515];
    
elseif(dz == 5)
    
    a = [9.51;11.17;12.48;13.6;14.3;13.7;13.5;15.8;14.7;11.8;11.4;10.9;10.5;10.1;9.7;9.3;8.96;8.6;8.3;7.95];
    b = -[25.9;107.5;232;344.6;457;433;520;1112;1125;853;931;1007;1080;1149;1216;1281;1343;1403;1460;1515];
    
    % if just duplicate, get steppy behavior
    % try to take means
    aa(1) = a(1);
    bb(1) = b(1);
    
    for i=1:length(a)-1
        
        aa(i+1)=nanmean(a(i:i+1));
        bb(i+1)=nanmean(b(i:i+1));
        
    end
    
    newa = nan*ones(40,1);
    newb = newa; 
    
    newa(1:2:end) = aa; 
    newa(2:2:end) = a;
    newb(1:2:end) = bb; 
    newb(2:2:end) = b;
    
    a = newa./2; 
    b = newb./2; 

end

as = 28.194; % for sail area
bs = -1420.2;    

% initialize arrays
icebergs.Z = dz:dz:500; icebergs.Z=icebergs.Z';
zlen = length(icebergs.Z);
temp = nan.*ones(zlen,length(L));  % 100 layers of 5-m each, so up to 500 m deep berg
temps = nan.*ones(1,length(L));  % sail area

% do icebergs K<200 first, using formula
ind = find(K<=200);
if(~isempty(ind))
    for i=1:length(ind)
       kz = K(ind(i)); %keel depth
       kza = ceil(kz./dz); % layer index for keel
       for nl=1:kza
          temp(nl,ind(i)) = a(nl).*L(ind(i)) + b(nl); 
       end
    end
    % now do sail area
    temps(ind) = as.*L(ind) + bs;
    temps(L<65) = 0.077.*L(L<65).^2; % fix for L<65, barker 2004
end

% then do icebergs D>200 for tabular
indt = find(K>200);
if(~isempty(indt))
    for i=1:length(indt)
       kz = K(indt(i)); %keel depth
       kza = ceil(kz./dz); % layer index for keel
       for nl=1:kza
          temp(nl,indt(i)) = L(indt(i)) .* dz; % 10 m times length 
       end
    end
    % now do sail area, using 0.9 for rho_i and 1.024 for rho_sw, then 12.1% of area(or volume, really) should be above water
    temps(indt) = 0.1211.* L(indt).*K(indt);
end

% now convert to output variables
icebergs.CA = temp;
icebergs.uwL = temp./dz; 

% now use L/W ratio of 1.62:1 (from Dowdeswell et al.) to get widths
icebergs.uwW = icebergs.uwL./LWratio; 

% now calculate volumes/mass of each layer, assuming rectangular
dznew = dz.*ones(size(icebergs.uwL));
% for i=1:length(K)
%    kind = round(K(i)./dz);
%    dz(kind,i) = abs(K(i)-kind.*10); % bottom layer is not dz thick
%    dz((kind+1):end,i) = nan;
% end
icebergs.uwV = dznew.*icebergs.uwL.*icebergs.uwW;  
%icebergs.uwM = 900.*icebergs.uwV;  % in kg

% before, I would calculate sail areas and volumes using formula, but this leads
% to icebergs not in the right ratio in terms of rho_i / rho_w