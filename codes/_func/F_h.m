function fh = F_h(K,H,ALFA)

% F_H.M 
% -------------------------------------------------------------------------
%   CM Cobb-Douglas production function function: marginal product of labor
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also F, F_K

fh = K^ALFA*H^(1-ALFA)*(1-ALFA)/H;