% SETTINGS.M
%--------------------------------------------------------------------------

global approx sigma

    approx = 1; % Order of Taylor series expansion. NOTE: max(approx) = 2.
    sigma = 1;  % Shock covariance scaling parameter for SGU algorithm
    
% % OPTIONS:
% %=========


SWEEP1 = 'off';       % Set to 'on' if running batch parameter sweep; 'off' 
                      % if just running this file.

                      
UNIX = 1;             % System type
DO_IMPULSE1 = 1;      % Compute impulse response from linear model
DO_PLOT = 1;
MULTIPLOT = 1;        % Separate IRF plots
SHOCKUNIT = 1;        % Choice of impulse shock: 1% (0.01 unit shock) == 1
                      %                                  1 s.d. shock == 0.
MOMENTS = 1;          % Calculate moments
DATA = 0;             % Availability of DATA statistics ( ==1).

GLOBAL_DISPLAY = 'on';

    % Setting current path and adding subpaths
    cwd_str = cd;
    
    if UNIX ~= 0 % for UNIX/Linux/Mac users
        addpath(genpath(strcat(cwd_str,'/','_func/')));
        addpath(genpath(strcat(cwd_str,'/','_mat/')));
        addpath(genpath(strcat(cwd_str,'/','_tools/')));  
        addpath(genpath(strcat(cwd_str,'/','_latex/')));
        addpath(genpath(strcat(cwd_str,'/','_func/')));
    else % for PC-Windows users
        addpath(genpath(strcat(cwd_str,'\','_func\')));
        addpath(genpath(strcat(cwd_str,'\','_mat\')));
        addpath(genpath(strcat(cwd_str,'\','_tools\')));  
        addpath(genpath(strcat(cwd_str,'/','_latex\')));
    end