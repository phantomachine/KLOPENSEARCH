function out = demand(Plevel,Pi,Y,EPSILON,SHARE)

% DEMAND.M
% -------------------------------------------------------------------------
%   Demand for intermediate good derived from CES aggregator of
%   intermediate goods.
%    
% THETA   : Home bias share
% EPSILON : 1/EPSILON - 1 is elasticity of substition btw x-y
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See also GG, GG_1, GG_2
    
    if EPSILON == Inf
        out = Y*(SHARE*Plevel/Pi);
    elseif EPSILON == 1
        warning('cpi:Error','Limit EPSILON = 1 not admissible!')
    elseif EPSILON == 0
        warning('cpi:Error','Leotief EPSILON = 0 not admissible!')
    elseif EPSILON ~= Inf || EPSILON ~= 0 || EPSILON ~= 1
        out = Y*(SHARE*Plevel/Pi)^(EPSILON/(EPSILON-1));
    end
    
    