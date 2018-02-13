% % PARETO_RUN.m
% % -----------------------------------------------------------------------
% % This SCRIPT executes the model.
% % 
% % INPUT: none
% % 
% % OUTPUT: model -- a CLASS containing OBJECTS: 
% %             .x          : current states
% %             .xp         : next-period states
% %             .y          : current control
% %             .yp         : next-period control
% %             .f          : vector of equilibrium equations at zero
% %             .anal_deriv : sub-CLASS containing symbolic partial 
% %                           derivatives of f(x,y,xp,yp)
% %             .ss         : contains steady state values
% %             .param      : contains parameters
% %
% % REQUIRES: MATLAB Symbolic Math Toolbox
% %
% % -----------------------------------------------------------------------
% % Acknowledgement:
% % Based on algorithms by (c) Stephanie Schmitt-Grohe and Martin Uribe
% % -----------------------------------------------------------------------
% % (c), 2009 - , Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------

%close all
%clear 
%clc

% % =======================================================================
% % STEP 0: Global Settings ...
% %         REMARK: "approx" is User-defined
% % =======================================================================
% close all
% clear all
% clc

RUN_OLD = 0;

global eta

settings;

    
% % =======================================================================
% % STEP 1: Parse model equilibrium and variable definitions ...
% %         REMARK: User-defined
% % =======================================================================
disp(' ');
disp('STEP 1-2: Interpreting model equilibrium conditions')
model = pareto_eq;

xp = model.xp;
yp = model.yp;
x = model.x;
y = model.y;
f = model.f;
                                                           
nx = size(x,2);
ny = size(y,2);

% % =======================================================================
% % STEP 2-3: Solve model.f for deterministic steady state variables ...
% %           NOTES: Example has no closed form to f(x,y) = 0. 
% %                REMARK: User-defined. (Need numerical solver.)
% % =======================================================================

    disp(' ');                                         
    disp('STEP 2-3: Solving now for nonstochastic steady state levels')                                       
        pareto_ss;
    disp('STEP 2-3: Done!')

    nshock = size(model.params{16},1);
    model.options.NSHOCK = nshock;


eta = zeros(nx,nshock);


OMEGA = model.params{16};              % Exogenous shock covariance matrix

eta(end-nshock+1:end,:) = OMEGA;      % nshock exogenous 
                                      % processes are ordered 
                                      % in the last nshock rows 
                                      % of the state vector
if RUN_OLD == 0
% % =======================================================================    
% % STEP 4: Find all analytical partial derivatives of f at [x,y,xp,yp] ...
% %         Output "model" contains a new sub-CLASS: model.anal_deriv
% % =======================================================================
    disp(' ');
    disp('STEP 4: Symbolic derivatives ... Please wait')
    model.f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));   % f with change of 
                                                        % variables, e.g.
                                                        % x := log(x), etc.

    model = anal_deriv(model,approx);
    disp('STEP 4: Done!')
    

% % =======================================================================
% % STEP 5: Use STEP 4 solutions to parameterize STEP 3 derivatives
% % =======================================================================
    
    disp(' ');                                         
    disp('STEP 5: Numerical evaluation of derivatives')
    
    % Parameters
           [ ALPHA      ... % 1
             BETA       ... % 2
             GAMA       ... % 3
             DELTA      ... % 4
             PHI        ... % 5
             ETA        ... % 6
             SIGMA      ... % 7
             THETA      ... % 8
             EPSILON    ... % 9
             ZBAR       ... % 10
             A          ... % 11
             B          ... % 12
             C          ... % 13 
             b          ... % 14
             RHO        ... % 15
             ECOV       ... % 16
             OMEGA_I    ... % 17
             OMEGA_F    ... % 18
             LAM_1      ... % 19
             LAM_2      ... % 20
             LAM_3      ... % 21
             LAM_4      ... % 22
             LAM_5      ... % 23
             CHI        ... % 24
             SG                     ] = model.params{:};  % 25

           [ k1        ... % 1
             X1        ... % 2
             K1        ... % 3
             H1        ... % 4
             q1        ... % 5
             , ~,        ... % 6
             yf1       ... % 7
             yh1       ... % 8
             G1        ... % 9
             Z1        ... % 10
             Y1        ... % 11
             I1        ... % 12
             RER1                    ] = model.steady1{:};  % 13

           [ k2        ... % 1
             X2        ... % 2
             K2        ... % 3
             H2        ... % 4
             q2        ... % 5
             , ~,        ... % 6
             yf2       ... % 7
             yh2       ... % 8
             G2        ... % 9
             Z2        ... % 10
             Y2        ... % 11
             I2        ... % 12
             RER2                    ] = model.steady2{:};   % 13

    % Unpacking RHO
    
    RHO_Z1 = RHO(1,1);
    RHO_Z12 = RHO(1,2);
    RHO_Z21 = RHO(2,1);
    RHO_Z2 = RHO(2,2);
    
    % Log Steady state
    K1 = log(K1); 
    K2 = log(K2);
    Z1 = log(Z1);
    Z2 = log(Z2);
    
    X1 = log(X1);
    X2 = log(X2);
    H1 = log(H1);
    H2 = log(H2);
    q1 = log(q1);
    q2 = log(q2);
    yh1 = log(yh1);
    yf1 = log(yf1);
    yf2 = log(yf2);
    yh2 = log(yh2);
    Y1 = log(Y1);
    Y2 = log(Y2);
    RER1 = log(RER1);
    I1 = log(I1);
    I2 = log(I2);
    TOTX = log(TOTX);
    RELX = log(RELX);
    YPI1 = log(YPI1);
    YPI2 = log(YPI2);
    
    % Deflated by GDP deflator: 15/8/2010
    GDP1 = log(GDP1);
    GDP2 = log(GDP2);
    CON1 = log(CON1);
    CON2 = log(CON2);
    INV1 = log(INV1);
    INV2 = log(INV2);
    
    K1p = K1; 
    K2p = K2;
    Z1p = Z1;
    Z2p = Z2;
    
    X1p = X1;
    X2p = X2;
    H1p = H1;
    H2p = H2;
    q1p = q1;
    q2p = q2;
    yh1p = yh1;
    yf1p = yf1;
    yf2p = yf2;
    yh2p = yh2;
    Y1p = Y1;
    Y2p = Y2;
    TOTXp = TOTX;
    RELXp = RELX;
    RER = RER1;
    RERp = RER;
    I1p = I1;
    I2p = I2;
    YPI1p = YPI1;
    YPI2p = YPI2;

    
    % Un-pack output from model.anal_deriv sub-CLASS
      fx = model.anal_deriv.fx;
      fxp = model.anal_deriv.fxp;
      fy = model.anal_deriv.fy;
      fyp = model.anal_deriv.fyp;
      fypyp = model.anal_deriv.fypyp;
      fypy = model.anal_deriv.fypy;
      fypxp = model.anal_deriv.fypxp;
      fypx = model.anal_deriv.fypx;
      fyyp = model.anal_deriv.fyyp;
      fyy = model.anal_deriv.fyy;
      fyxp = model.anal_deriv.fyxp;
      fyx = model.anal_deriv.fyx;
      fxpyp = model.anal_deriv.fxpyp;
      fxpy = model.anal_deriv.fxpy;
      fxpxp = model.anal_deriv.fxpxp;
      fxpx = model.anal_deriv.fxpx;
      fxyp = model.anal_deriv.fxyp;
      fxy = model.anal_deriv.fxy;
      fxxp = model.anal_deriv.fxxp;
      fxx = model.anal_deriv.fxx;
 
    % Evaluate above numerically
    num_eval;
    
% % =======================================================================
% % STEP 6: Solve for stable REE dynamics by SGU pertubation method.
% % =======================================================================

    disp(' ');                                         
    disp(['STEP 6: Solving for REE. Policy function approx. order = ',...
                                                          int2str(approx)]) 
                                                      
    if approx == 1
        %First-order approximation
        [gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
        
        if PHI > 1.01
            save pareto_solve_k gx hx
        else
            save pareto_solve gx hx
        end
        
    elseif approx == 2
        %First-order approximation
        [gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
        
               
        %Second-order approximation
        [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,...
                        nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,...
                        nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 

        [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,...
                        nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,...
                        nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
        if PHI > 1.01            
            save pareto_solve_k gx hx gxx hxx
        else
            save pareto_solve gx hx gxx hxx
        end
    end
    
    disp(' ') 
    disp('Check for stability')
    disp('eig(hx) = ')
    disp(eig(hx))
    disp(' ')
    
    disp('STEP 6: Done!')

elseif RUN_OLD == 1
    
    load pareto_solve
    
end
% % =======================================================================
% % STEP 7: Analyze dynamics
% % =======================================================================

                  VARNAMES = 'X';                            % 1                               
VARNAMES = char(VARNAMES, 'X^{\ast}');                       % 2
VARNAMES = char(VARNAMES, 'H');                              % 3
VARNAMES = char(VARNAMES, 'H^{\ast}');                       % 4
VARNAMES = char(VARNAMES, 'q');                              % 5
VARNAMES = char(VARNAMES, 'q^{\ast}');                       % 6
VARNAMES = char(VARNAMES, 'y_{h}');                          % 7
VARNAMES = char(VARNAMES, 'y_{f}');                          % 8
VARNAMES = char(VARNAMES, 'y_{f}^{\ast}');                   % 9
VARNAMES = char(VARNAMES, 'y_{h}^{\ast}');                   % 10
VARNAMES = char(VARNAMES, 'Y');                              % 11
VARNAMES = char(VARNAMES, 'Y^{\ast}');                       % 12
VARNAMES = char(VARNAMES, 'I');                              % 13
VARNAMES = char(VARNAMES, 'I^{\ast}');                       % 14
VARNAMES = char(VARNAMES, 'TOTX');                           % 15
VARNAMES = char(VARNAMES, 'RER');                            % 16
VARNAMES = char(VARNAMES, 'CON/CON^{\ast}');                 % 17
VARNAMES = char(VARNAMES, 'GDP');                            % 18
VARNAMES = char(VARNAMES, 'GDP^{\ast}');                     % 19
VARNAMES = char(VARNAMES, 'CON');                            % 20
VARNAMES = char(VARNAMES, 'CON^{\ast}');                     % 21
VARNAMES = char(VARNAMES, 'INV');                            % 22
VARNAMES = char(VARNAMES, 'INV^{\ast}');                     % 23
VARNAMES = char(VARNAMES, 'P_{Y}');                          % 24
VARNAMES = char(VARNAMES, 'P_{Y}^{\ast}');                   % 25

VARNAMES = char(VARNAMES, 'K');                              % 26
VARNAMES = char(VARNAMES, 'K^{\ast}');                       % 27
VARNAMES = char(VARNAMES, 'TFP');                            % 28
VARNAMES = char(VARNAMES, 'TFP^{\ast}');                     % 29


% % 7(a) Impulse response functions:
% % -----------------------------------------------------------------------
model.out.VARNAMES = VARNAMES;
model.out.gx = gx;
model.out.hx = hx;

% % Options:
model.options.IMPSELECT = [24,25,16,26,22,20,5,6];
model.options.SHOCKSELECT = 1:2;
model.options.HORIZON = 30;
model.options.LINETYPE = {'+k'}; %{'or'};

if DO_IMPULSE1 == 1 && DO_PLOT == 1
    model.options.DO_PLOT = 1;
    
else
    model.options.DO_PLOT = 0; 
end

    if MULTIPLOT == 1
        model.options.MULTIPLOT = 1;
    else
        model.options.MULTIPLOT = 0;
    end

    if SHOCKUNIT == 1
        SHOCKCOV = 0.01*eye(size(OMEGA)); % 1% shock
    else
        SHOCKCOV = sqrt(OMEGA);            % 1 std. dev. shock
    end

model.options.x_init_mat = [zeros(nshock,nx-nshock), SHOCKCOV ];   

    if DO_IMPULSE1 == 1
        disp(' ')
        disp('STEP 7(a): Compute and plot Impulse responses')
        disp(' ')
        [ IRMAT ] = impulse_linear(model);
    end
    
% % % 7(b) Autocovariograms:
% % % -----------------------------------------------------------------------

if MOMENTS == 1
        disp(' ')
        disp('STEP 7(b): Compute moments')
        disp(' ')
        
    varshock = blkdiag(zeros(nx-nshock,nx-nshock), OMEGA);
    
    % % Options:
    model.options.method = 0; % 0 == algebraic inversion; 1 == doubling algorithm
    model.options.J = 1;
    
    model.options.MOM_SELECT_y = 1:23;                   %  1 : 14
    model.options.MOM_SELECT_x = [];                       % 17 : 22
    
    model.options.Y_INDEX = 18;
    model.options.NX1_INDEX = [];
    model.options.NX2_INDEX = [];
    
    % Cross-country correlation selector matrix: corr(x, x*)
    model.options.CROSSCOUNTRY = [  18, 19;                % (Y,Y*)
                                    20,  21;               % (CON,CON*)
                                    16, 17];               % (RER,X/X*)
        
    % Generate moments tables:
    [STD_table,XCF_Y_table,modelcorrcc,modelautocorr] = ...
                            mom_out(model,DATA,varshock,GLOBAL_DISPLAY);

end

disp('--------- PARETO_RUN: THE END ----------------------------')
    