% PARAMETERSCOMMON.M
 
%-------------------------------------------------------------------------
% 1. COMMON PARAMETERS TO ALL MODELS

% MICRO PARAMETERS
%ALPHA = 0.288;   % Income share of Capital in F(K,H)
BETA = 0.99;      % Discount factor
GAMA = 1;         % CRRA parameter: CM good X
DELTA = 0.1/4;    % Capital depreciation rate (Heathcote-Perri)
%PHI = 1.2;       % Capital complementarity in DM technology: c(q,K)
ETA = 1.0;        % CRRA parameter: DM good q
%SIGMA = 0.26/4;  % Pr{Match} = Pr{OR(buyer,seller)}, (AWW's)/4
THETA = 0.9397;   % Home bias in intermediate good
KAPPA =  0.85;    % Proportion of money transactions in DM: Cooley-Hansen (1991)
 
% Calibration Targets
KY_target = 2.23*4;%3.23; %3.23;
v_target = 5.29/4; % M1 velocity (quarterly average)
H_target = 1/3;
LS_target = 0.712;
NTS_target = 0.5;
MUP_target = 0.3;
MDE_target = -0.23;
MKP_target = 0.33;

% Elasticity of subs H-F, EPSILON/(1-EPSILON) = 1.5: BKK (1995)
SIG_epsilon = 1.5; 
        EPSILON = SIG_epsilon/(SIG_epsilon-1);
        
ZBAR = 1;        % Steady state TFP
EPSIBAR = 1;     % Steady state gross money supply growth factor

%A = 13;         % Disutility of labor: CM
%B = 1;          % Scale parameter on CM CRRA utility
C = 1;           % Scale parameter on DM CRRA utility
b = 0.0001;      % Bounding from zero consumption parameter

SG = 0.0;        % Share of govt consumption in output G1, G2

% COMPOSITE PARAMETERS
OMEGA_I = (THETA + (1-THETA)*(THETA/(1-THETA))^(1/(1-EPSILON)))^EPSILON;

OMEGA_F = 1 + (THETA/(1-THETA))^(EPSILON/(1-EPSILON));

% STOCHASTIC PROCESSES - TFP

% % VAR(1) - Heathcote-Perri Estimates
% RHO =[ 0.92     0.025     ;
%       0.025       0.92   ]; 

    % VAR(1) - CKM Calibration
    RHO =[ 0.95     0     ;
           0       0.95   ];

    % Variance-Covariance Matrix of TFP shocks - CKM
    ECOV = (0.007^2)*[  1       0.25    ;
                        0.25    1       ]; 
                    
    % % Shlagenhauf-Wrase (JIMF, 1995)       
    % ECOV = [ 0.01359     0.00191     ;
    %          0.00191     0.01261 ]; 
 
%-------------------------------------------------------------------------
% 2. COMMON PARAMETERS TO {COMPECON, BONDECON}
%-------------------------------------------------------------------------

% STOCHASTIC PROCESSES - Money growth
mgrowth_process = 'swr';  % { 'ckm', 'swr', 'data' }

if strcmp(mgrowth_process,'ckm') == 1
    disp('PARAMETERCOMMON.M:')
    disp(' ')
    disp('Chosen Money Supply process: Chari-Kehoe-McGrattan (ReStud, 2003)')
    
    MU = [ 0.68    0.0;
          0.00    0.68 ];

    MCOV = [ 0.023    0.000    ;
             0.000    0.023   ].^2;
         
elseif strcmp(mgrowth_process,'swr') == 1 % Shlagenhauf-Wrase (JIMF, 1995)
    disp('PARAMETERCOMMON.M:')
    disp(' ')
    disp('Chosen Money Supply process: Shlagenhauf-Wrase (JIMF, 1995)')
    
    MU = [ 0.592    0.007;
          0.098    0.695 ];
    
    MCOV = [ 0.00397    0.00026    ;
             0.00026    0.00662   ];
elseif strcmp(mgrowth_process,'data') == 1
    disp('PARAMETERCOMMON.M:')
    disp(' ')
    disp('Chosen Money Supply process: Data estimate, VAR(1)')
    disp('Option currently not available')
    break
end


     
% Tax rates
TAU_X = 0.069;
TAU_K = 0.548;
TAU_H = 0.242;