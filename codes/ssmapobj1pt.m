function f = ssmapobj1pt(x,~,~,~,~,~,~,...
                        ~,~,~,~,~,~,~, ...
                              ~, ~, ~, ...
                                H_target,KY_target,NTS_target) 

% --------------------------------------------------------------------
% SSMAPOBJ.M
% --------------------------------------------------------------------
% Quadratic Criterion Function for calibration targets. 
% --------------------------------------------------------------------
% (c) 2010- T.Kam; Email: mortheus@gmail.com    
%
% See also SSMAPCON

% Extract relevant candidate variables:                            
    H = x(5);        % in (0,1)
    KYratio = x(6);  % > 0
    NTS = x(7);     % in (0,1)

% Quadratic loss function on target steady-state mean variables:
f = (H/H_target -1)^2 ...                    % Labor
        + (KYratio/KY_target - 1)^2 ...       % Capital-output ratio
                        + (NTS/NTS_target -1)^2;     % Nontradable share in 
                                                  % total consumption