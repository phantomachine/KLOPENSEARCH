function out = GG(x,y,EPSILON,THETA)

% GG.M
% -------------------------------------------------------------------------
% CES aggregator of intermediate goods.
%    
% THETA   : Home bias share
% EPSILON : 1/EPSILON - 1 is elasticity of substition btw x-y
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See also GG_2, GG_1, DEMAND

    if EPSILON == Inf
        out = (x^THETA)*(y^(1-THETA));
    elseif EPSILON == 1
        out = x*THETA + y*(1-THETA);
    elseif EPSILON == 0
        warning('GG:Error','Leotief EPSILON = 0 not admissible!')
    elseif EPSILON ~= Inf || EPSILON ~= 0 || EPSILON ~= 1
        out = (THETA*x^(1/EPSILON) + (1-THETA)*(y^(1/EPSILON)))^EPSILON;
    end