function [IMPMAT] = imprespfun(A,C,G,HORIZON,SHOCKUNIT,OMEGA)

% % Copyleft 2006, T.Kam
% Linear state space system: x' = Ax + Cv' and y = Gx.
%
% Inputs:
%
% Matrices: A, C, G
% HORIZON: Length of impulse response functions
% SHOCKUNIT: ~0 == create SHOCKUNIT s.d. of shock
%            =0 == create 1 unit of shock
% OMEGA: Covariance matrix of shocks, E(vv')

nx = size(A,1);
ny = size(G,1);
nshock = size(C,2); % number of primitive iid shocks

IMPMAT = zeros(nx+ny,HORIZON+1,nshock); % (nx + ny) x T x nz dimensional array

for j = 1:nshock
SHOCK = zeros(1,nshock);    
    
if SHOCKUNIT == 0
        SHOCK(j) = 1;
    else
        SHOCK(j) = SHOCKUNIT*sqrt(OMEGA(j,j));  % SHOCKUNIT x 1 s.d. shock
end

    x = C*SHOCK';
    y = G*x;
    
    for i = 1:HORIZON
        xp = A*x(:,i);
        y = [ y, G*xp];
        x = [ x, xp ];
    end   
    
IMPMAT(:,:,j) = [x; y];
end
    


