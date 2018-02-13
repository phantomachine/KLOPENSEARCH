function [ g_kval ] = g_k(q, K, Z, XI)

% G.M
% -------------------------------------------------------------------------
% g_k derivative function in price taking version.
%
% -------------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See G, G_Q, COST, C_Q

g_kval = (XI/(1+XI)) * (q/Z) * c_qk(q/Z,K,XI) /Z;