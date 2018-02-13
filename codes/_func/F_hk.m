function fhk = F_hk(K,H,ALFA)

% F_HK.M 
% -------------------------------------------------------------------------
%   CM Cobb-Douglas production function function: Cross-partial (H,K)
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also F, F_K, F_H, F_HH

fhk = ALFA*(1-ALFA)*(K^(ALFA-1))*H^(-ALFA);