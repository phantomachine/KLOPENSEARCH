function [ gnash_q_val ] = g_nash_q(q,K,Z,BTHETA,PHI,ETA,b,C)

% G_NASH_Q.M
% -------------------------------------------------------------------------
% Derivative function g_q in proportional bargaining version.
% BTHETA is bargaining power of buyer.
% -----------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
%
% See also U_Q, C_Q, U_QQ, C_QQ, G_NASH, GAMMA_NASH

% Primitive value
u = uq(q,ETA,b,C);
c = cost(q/Z,K,PHI);

% First derivative value
uqd = u_q(q,ETA,b,C);
cq = c_q(q/Z,K,PHI)/Z;

% Second derivative value
uqq = u_qq(q,ETA,b,C);
cqq = c_qq(q/Z,K,PHI)/(Z^2);

% Construct g_q(q,K,Z) 

numerator = ( uqd * cq )*(BTHETA*uqd + (1-BTHETA)*cq  ) ...
                + BTHETA*(1-BTHETA)*(u - c)*(uqd*cqq - cq*uqq);
                                
denominator = ( BTHETA*uqd + (1-BTHETA)*cq )^2;

gnash_q_val = numerator / denominator;