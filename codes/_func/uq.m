function [uqval] = uq(q, ETA2, b, KONSTANT)

% UQ.M 
% -------------------------------------------------------------------------
%   CRRA utility of q
% -----------------------------------------------------------------------
%   (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% -----------------------------------------------------------------------
% 
%See also U_Q
    
%if ETA2 > 1
%       uqval = KONSTANT*(( q + b ).^(1 - ETA2) - b^(1 - ETA2))/(1 - ETA2);
%elseif ETA2 == 1
        uqval = KONSTANT*( log(q + b) - log(b) );
%end