function fkk = F_kk(K,H,ALFA)

% F_H.M 
% -------------------------------------------------------------------------
%   CM Cobb-Douglas production function function: second derivative wrt K
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also F, F_H

fkk = ALFA*(ALFA-1)*K^(ALFA-2)*H^(1-ALFA);