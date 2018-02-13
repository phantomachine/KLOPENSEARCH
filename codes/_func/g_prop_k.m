function [ g_p_kval ] = g_prop_k(q,K,Z,BTHETA,PHI)

% G_PROP_K.M
% -------------------------------------------------------------------------
% Derivative function g_k in proportional bargaining version.
% BTHETA is bargaining power of buyer.
% -----------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
%
% See also U_Q, C_Q, G_PROP, GAMMA_PROP


g_p_kval = BTHETA*c_k(q/Z,K,PHI);