function K = keeldepth(L,method)
%
% calculate keel depth based on waterline length based on empirical formula
%
% L is vector of waterline lengths in meters, round to nearest 10
%
% methods are 'auto','barker','hotzel','constant'
%   (auto uses length to switch, while constant is K=0.7L)
% for icebergs with L<=160, use Barker et al. 2004 formula, K= 2.91L^0.71
% 
% for icebergs with L>=170, assume they are tabular, use Hotzel&Miller, K = 3.781L^0.63
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% round L to nearest 10 m, to decide which formula to use
a = round(L./10);
La = 10.*a;

bind = La<=160; % use barker
hind = La>160;  % use hotzel

if nargin==1
    method = 'auto';
end

switch method
    case 'auto'  % use both
        K(hind) = 3.78 .* La(hind).^0.63;
        %
        K(bind) = 2.91 .* La(bind).^0.71; 
    
    case 'barker'  % use only barker
        K = 2.91 .* La.^0.71; 
    
    case 'hotzel'  % use only hotzel
        K = 3.78 .* La.^0.63;
        
    case 'constant'  % use K = 0.7L
        K = 0.7.*La;
    
    case 'mean'  % do all methods and take mean
        temp = nan.*ones(length(L),4);
        temp(bind,1) = 2.91 .* La(bind).^0.71; 
        temp(hind,1) = 3.78 .* La(hind).^0.63;
        temp(:,2) = 2.91 .* La(:).^0.71;
        temp(:,3) = 3.78 .* La(:).^0.63;
        temp(:,4) = 0.5.*La(:);
        % take mean
        K = mean(temp,2); 
end
