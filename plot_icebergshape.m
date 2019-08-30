function [axm,shape] = plot_icebergshape(L, uwL, freeB, Z, K, col) 
% plot iceberg shape based on variables in ice structure,
%  like from init_iceberg_size.m mfile
%
% input: L = length
%            uwL = vector of underwater lengths
%            freeB = freeboard
%            Z = vector of layer depths
%            K = keel depth
%            col = color to use
%
% output   axm = maximum of L/2,K,freeB,max UWL/2
%          shape = polygon to plot shape of iceberg without layers
%                  where shape(:,1) is x, shape(:,2) is y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz = diff(Z(1:2));   % thickness of layers
keeli = ceil(K./dz); % index of keel
dzk = -((keeli-1).*dz - K); %thickness of last layer

% plot above water portion
xx = [-L/2 L/2 L/2 -L/2 -L/2];
yy = [0 0 freeB freeB 0];

xk = min(xx).*ones(2,1); %shape of iceberg
yk = [max(yy);min(yy)]; 

plot(xx,yy,'linewidth',1.25,'color',col);
   
hold on; 
grid on;

% now plot layers   
for i = 1:keeli-1
    
  dt = uwL(i) ./ 2;
  
  xx = [-dt dt dt -dt -dt];
  yy = -[Z(i)-dz Z(i)-dz Z(i) Z(i) Z(i)-dz];  

  xk = [xk; min(xx) .* ones(2,1)]; %shape of iceberg
  yk = [yk; [max(yy); min(yy)]];
  
  plot(xx,yy,'linewidth',1.25,'color',col);

end

% now do bottom layer
dt = uwL(keeli)./2;
xx = [-dt dt dt -dt -dt];
yy = -[Z(keeli)-dz Z(keeli)-dz Z(keeli-1)+dzk Z(keeli-1)+dzk Z(keeli)-dz];

xk = [xk;min(xx).*ones(2,1)]; %shape of iceberg
yk = [yk;[max(yy);min(yy)]];
 
shape(:,1) = [xk; flipud(-xk); xk(1)];
shape(:,2) = [yk; flipud(yk); yk(1)];

plot(xx,yy,'linewidth',1.5,'color',col);

axm = [L/2  K  freeB  nanmax(uwL(:))./2];  
