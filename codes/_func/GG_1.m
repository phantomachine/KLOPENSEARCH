function out = GG_1(x,y,EPSILON,THETA)
    
% GG_1.M
% -------------------------------------------------------------------------
%   Partial Derivative of CES aggregator of intermediate goods 
%   wrt x (1st argument).
%    
% THETA   : Home bias share
% EPSILON : 1/EPSILON - 1 is elasticity of substition btw x-y
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See also GG, GG_2, DEMAND
    
    if EPSILON == Inf
        out = THETA*(x/y)^(THETA-1);  % Cobb-Douglas MP_x
    elseif EPSILON == 1
        out = THETA;                  % Pefect substitutes, constant MP_x
    elseif EPSILON == 0
        warning('GG:Error','Leotief EPSILON = 0 not admissible!')
    elseif EPSILON ~= Inf || EPSILON ~= 0 || EPSILON ~= 1
        out = THETA*x^(1/EPSILON - 1)...
          *(THETA*x^(1/EPSILON) + (1-THETA)*(y^(1/EPSILON)))^(EPSILON - 1);
                                      % CES case for MP_x
    end