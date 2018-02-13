function f = ssmapobj3propbarg(x,~,~,~,~,...
                        ~,~,~,~,~,~,~, ...
                              ~, ~, ~,~, ...
                                H_target, KY_target, v_target, ...
                                LS_target, NTS_target,MUP_target) 

% --------------------------------------------------------------------
% SSMAPOBJ.M
% --------------------------------------------------------------------
% Quadratic Criterion Function for calibration targets. 
% --------------------------------------------------------------------
% (c) 2010- T.Kam; Email: mortheus@gmail.com    
%
% See also SSMAPCON

% Extract relevant candidate variables:                            
    H = x(8);        % in (0,1)
    KYratio = x(9);  % > 0
    velocity = x(10); % > 0
    LS = x(11);       % in (0,1)
    NTS = x(12);     % in (0,1)
    MUP = x(13);     % > 0

% Quadratic loss function on target steady-state mean variables:
f = (H - H_target)^2 ...                    % Labor
        + (KYratio - KY_target)^2 ...       % Capital-output ratio
                + (velocity - v_target)^2 ...  % Velocity of Money
                    + (LS - LS_target)^2 ...        % Labor share
                        + (NTS - NTS_target)^2 ...  % Nontradable share in 
                         + (MUP - MUP_target)^2;        % total consumption