function cqk = c_qk(q,K,XI)

% C_QK.M 
% -------------------------------------------------------------------------
%   Cross-partial derivative of utility effort cost of q in the DM. Assume:
%
%   q = q(cost, K) = cost^(1/XI) * K^((XI-1)/XI),
%
% is Cobb-Douglas, where the share of K is (XI - 1)/XI and XI >= 1.
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also COST, C_K, C_Q, C_QQ

    cqk = (1-XI)*XI*( q^(XI-1) )*K^(-XI);