function [ gnash ] = g_nash(q,K,Z,BTHETA,PHI,ETA,b,C)

% G_NASH.M
% -------------------------------------------------------------------------
% g function in Nash Bargaining version. 
%
% -------------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See also GAMMA_PROP, G_PROP_Q, UTIL, U_Q, COST, C_Q

% Primitive value
u = uq(q,ETA,b,C);
c = cost(q/Z,K,PHI);

% First derivative value
uqd = u_q(q,ETA,b,C);
cq = c_q(q/Z,K,PHI)/Z;

% Construct g function
numerator = BTHETA * c * uqd + (1-BTHETA) * u * cq;
                                
denominator = BTHETA * uqd + (1-BTHETA) * cq;

gnash = numerator / denominator;