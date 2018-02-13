function [ gval ] = g_prop_q(q,K,Z,THETA_B,XI,ETA,b,C)

% G_PROP_Q.M
% -------------------------------------------------------------------------
% Derivative function g_q in proportional bargaining version.
% BTHETA is bargaining power of buyer.
% -----------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
%
% See also U_Q, C_Q, G_PROP, GAMMA_PROP


% gval = (1-THETA_B)*u_q(q,ETA,b,C) + THETA_B*c_q(q/Z,K,PHI)/Z;

gval = (1-THETA_B)*(C*(q+b)^(1-ETA)/(q+b)) + THETA_B*(XI*( (q/Z) / K )^(1-XI))/Z;
