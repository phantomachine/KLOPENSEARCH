function cq = c_q(q,K,XI)

% C_K.M 
% -------------------------------------------------------------------------
%   Marginal utility of effort cost of q in the DM. Assume:
%
%   q = q(cost, K) = cost^(1/XI) * K^((XI-1)/XI),
%
% is Cobb-Douglas, where the share of K is (XI - 1)/XI and XI >= 1.
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also COST, C_K

    cq = XI*( q / K )^(XI-1);