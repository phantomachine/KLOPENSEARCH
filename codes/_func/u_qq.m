function uqq = u_qq(q, ETA2, b, C)

% U_QQ.M 
% -------------------------------------------------------------------------
%   CRRA second derivative function for utility of q
% -----------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
% 
%See also UQ, U_Q

    uqq = -ETA2*C*(q+b)^(-1-ETA2);