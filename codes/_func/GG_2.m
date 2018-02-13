function out = GG_2(x,y,EPSILON,THETA)

% GG_2.M
% -------------------------------------------------------------------------
%   Partial Derivative of CES aggregator of intermediate goods 
%   wrt y (2nd argument).
%    
% THETA   : Home bias share
% EPSILON : 1/EPSILON - 1 is elasticity of substition btw x-y
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See also GG, GG_1, DEMAND

    if EPSILON == Inf
        out = (1-THETA)*(x/y)^(THETA);  % Cobb-Douglas MP_y
    elseif EPSILON == 1
        out = 1-THETA;                  % Perfect substitutes, const MP_y
    elseif EPSILON == 0
        warning('GG:Error','Leotief EPSILON = 0 not admissible!')
    elseif EPSILON ~= Inf || EPSILON ~= 0 || EPSILON ~= 1
        out = (1-THETA)*(y^(1/EPSILON - 1))...
          *(THETA*x^(1/EPSILON) + (1-THETA)*(y^(1/EPSILON)))^(EPSILON - 1);
                                        % CES MP_y
    end