function f = ssmapobj2pb(x,~,~,~,~,~,...
                        ~,~,~,~,~,~,~, ...
                              ~, ~,~,~,~, ...
                                H_target,KY_target,...
                                NTS_target,MKP_target) 

% --------------------------------------------------------------------
% SSMAPOBJ2.M
% --------------------------------------------------------------------
% Quadratic Criterion Function for calibration targets. 
% --------------------------------------------------------------------
% (c) 2010- T.Kam; Email: mortheus@gmail.com    
%
% See also SSMAPCON

%disp(MKP_target)

% Extract relevant candidate variables:                            
    
    H = x(6);        % in (0,1)
    KYratio = x(7);  % > 0
    NTS = x(8);     % in (0,1)
    MKP = x(9);     % in (0,m), m < +infty

% Quadratic loss function on target steady-state mean variables:
f = (H/H_target -1)^2                        ... % Labor
    + (KYratio/KY_target - 1)^2              ... % Capital-output ratio
    + (NTS/NTS_target -1)^2             ... % Nontradable consumption share
    + (MKP/MKP_target - 1)^2;               % Total markup