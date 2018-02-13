function cqq = c_qq(q,K,XI)

% C_QQ.M 
% -------------------------------------------------------------------------
%   Second derivative of utility effort cost of q in the DM. Assume:
%
%   q = q(cost, K) = cost^(1/XI) * K^((XI-1)/XI),
%
% is Cobb-Douglas, where the share of K is (XI - 1)/XI and XI >= 1.
% -------------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -------------------------------------------------------------------------
% 
%See also COST, C_K

    cqq = (XI-1)*XI*( q^(XI-2) )*K^(1-XI);