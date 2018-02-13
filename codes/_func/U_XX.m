function UX = U_XX(X,B,GAMMA)

% U_XX.M 
% -------------------------------------------------------------------------
%   CRRA utility of X: second derivative
% -----------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
% 
% See also U

    UX = -GAMMA*B*X^(-GAMMA-1);