function [ gammaval ] = gamma_nash(q,K,Z,BTHETA,PHI,ETA,b,C)

% GAMMA_NASH.M
% -------------------------------------------------------------------------
% "GAMMA" function in Nash bargaining version.
% BTHETA is bargaining power of buyer.
% -----------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
%
% See also G_PROP, G_PROP_Q, U_Q, C_K

% First derivative values
cq = c_q(q/Z,K,PHI)/Z;
ck = c_k(q/Z,K,PHI);

gk = g_nash_k(q,K,Z,BTHETA,PHI,ETA,b,C);
gq = g_nash_q(q,K,Z,BTHETA,PHI,ETA,b,C);

gammaval = ck - cq*gk/gq;