function fhh = F_hh(K,H,ALFA)

% F_HH.M 
% -------------------------------------------------------------------------
%   CM Cobb-Douglas production function function: second derivative wrt H
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also F, F_K

fhh = -ALFA*(1-ALFA)*(K^ALFA)*H^(-1-ALFA);