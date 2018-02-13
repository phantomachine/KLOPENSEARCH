function g_qval = g_q(q,K,Z,XI)

% G_Q.M
% -------------------------------------------------------------------------
% Derivative g_q function in price taking version.
%
% -------------------------------------------------------------------------
% (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
%
% See G, COST, C_Q

    g_qval = XI*c_q(q/Z,K,XI)/Z;