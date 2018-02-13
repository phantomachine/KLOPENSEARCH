 function [UX] = U(x, B, GAMMA)

% U.M 
% -------------------------------------------------------------------------
%   CRRA utility of X. Note GAMMA ~= 1
% -----------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
% 
% See also U_X

%    if ETA1 ~= 1
        UX = B * ( x^(1 - GAMMA) - 1 )/(1 - GAMMA) ;
%     elseif ETA1 == 1
%         UX = B * log(x) ;
%     end
end