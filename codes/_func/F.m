function [Fval] = F(K, H, ALFA)

% F.M 
% -------------------------------------------------------------------------
%   Provides CM Cobb-Douglas production function functions
% -----------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
% 
%See also F_H, F_K
    
Fval = ( K^ALFA ) * H^(1-ALFA) ;