function fk = F_k(K,H,ALFA)

% F_H.M 
% -------------------------------------------------------------------------
%   CM Cobb-Douglas production function function: marginal product of
%   capital.
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also F, F_H

fk = K^ALFA*ALFA/K*H^(1-ALFA);