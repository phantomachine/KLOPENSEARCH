% % PB_CALIBRA_RUN.m
% % -----------------------------------------------------------------------
% % This SCRIPT executes the Proportional Bargaining model with credit
% % and money in DM. It studies a deviation indexed by THETA_B from the
% % baseline Price Taking calibration, and what this deviation does
% % relative to the Price Taking model dynamics.
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

close all
clear 
clc


% % =======================================================================
% % STEP 0: Global Settings ...
% %         REMARK: "approx" is User-defined
% % ======================================================================= 

RUN_OLD = 1;

protocol = 'pb';

global eta

settings;

    
% % =======================================================================
% % STEP 1: Parse model equilibrium and variable definitions ...
% %         REMARK: User-defined
% % =======================================================================
disp(' ');
disp('STEP 1-2: Interpreting model equilibrium conditions')
model = pb_eq;

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
        pb_calibra_ss;
    disp('STEP 2-3: Done!')
    
    % dim(ECOV) + dim(MCOV)... See bondecon_ss.m :
    nshock = 4;
    model.options.NSHOCK = nshock;

eta = zeros(nx,nshock);

% Exogenous shock covariance matrix:
OMEGA = blkdiag(ECOV,MCOV);  
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
    
% Convert steady state values to log
    Z1 = log(Z1);    
    Z2 = log(Z2);
    
    psi1 = log(psi1);
    psi2 = log(psi2);
    
    K1 = log(K1);    
    K2 = log(K2);
       
    X1  = log(X1);   
    X2 = log(X2);
    
    H1 = log(H1);    
    H2 = log(H2);
    
    q1  = log(q1);   
    q2 = log(q2);
    
    qz1 = log(qz1);
    qz2 = log(qz2);
    
    yh1 = log(yh1);   
    yf2 = log(yf2);
    
    yf1 = log(yf1);   
    yh2 = log(yh2);
    
    %NX1 = log(NX1); Comment: Leave NX1 and NX2 in levels. 
    %NX2 = log(NX2);          Possible nonpositive realizations
    
    Y1 = log(Y1);   
    Y2 = log(Y2);
    
    I1 = log(I1);   
    I2 = log(I2);
    
    P1 = log(P1);   
    P2 = log(P2);
    
    Pyh1 = log(Pyh1);
    Pyf2 = log(Pyf2);
    
    Pyf1 = log(Pyf1);
    Pyh2 = log(Pyh2);
    
    RER = log(RER);
    E = log(E); 
    RELX = log(RELX);
    
    Ptilde1 = log(Ptilde1); %added 27/05/2010
    Ptilde2 = log(Ptilde2);
    YPI1    = log(YPI1);
    YPI2    = log(YPI2);
    RER_YPI = log(RER_YPI);

    Agg_C1  = log(Agg_C1);
    Agg_C2  = log(Agg_C2);

    Agg_Y1  = log(Agg_Y1);
    Agg_Y2  = log(Agg_Y2);

    Agg_I1  = log(Agg_I1);
    Agg_I2  = log(Agg_I2);
    
    P2P1 = log(P2P1);

    % Unpacking RHO -  TFP VAR(1)
    
    RHO_Z1 = RHO(1,1);
    RHO_Z12 = RHO(1,2);
    RHO_Z21 = RHO(2,1);
    RHO_Z2 = RHO(2,2);
    
    % Unpacking MU - Money growth VAR(1)
    
    MU_M1 = MU(1,1);
    MU_M12 = MU(1,2);
    MU_M21 = MU(2,1);
    MU_M2 = MU(2,2);
    
    % Steady state
    K1p = K1; 
    K2p = K2;
    Z1p = Z1;
    Z2p = Z2;
    psi1p = psi1;
    psi2p = psi2;
    
    X1p = X1;
    X2p = X2;
    H1p = H1;
    H2p = H2;
    q1p = q1;
    qz1p = qz1;
    q2p = q2;
    qz2p = qz2;
    yh1p = yh1;
    yf1p = yf1;
    yf2p = yf2;
    yh2p = yh2;
    Y1p = Y1;
    Y2p = Y2;
    I1p = I1;
    I2p = I2;
    NX1p = NX1;
    NX2p = NX2;
    P1p = P1;
    P2p = P2;
    Pyh1p = Pyh1;
    Pyf2p = Pyf2;
    Pyf1p = Pyf1;
    Pyh2p = Pyh2;
    RERp = RER;
    Ep = E;
    RELXp = RELX;
    TOTp = TOT;
    Ptilde1p = Ptilde1;
    Ptilde2p = Ptilde2;
    YPI1p    = YPI1;
    YPI2p    = YPI2;
    RER_YPIp = RER_YPI;

    Agg_C1p  = Agg_C1;
    Agg_C2p  = Agg_C2;

    Agg_Y1p  = Agg_Y1;
    Agg_Y2p  = Agg_Y2;

    Agg_I1p  = Agg_I1;
    Agg_I2p  = Agg_I2;
    
    P2P1p = P2P1;
    
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
        
        save pb_calibra_solve gx hx
        
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
        
        save pb_calibra_solve gx hx gxx hxx
        
    end
    
    disp(' ') 
    disp('Check for stability')
    disp('eig(hx) = ')
    disp(eig(hx))
    disp(' ')
    
    disp('STEP 6: Done!')
elseif RUN_OLD == 1
   
   load pb_calibra_solve
   
end

% % =======================================================================
% % STEP 7: Analyze dynamics
% % =======================================================================
    model.out.VARNAMES = var_names_common;
    model.out.gx = gx;
    model.out.hx = hx;
    
% % 7(a) Impulse response functions:
% % -----------------------------------------------------------------------

    % % Options:
        model.options.HORIZON = 30;
        model.options.LINETYPE = {'-k'}; %{'or'};

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

    for figno = 1:2 % loop for 2 separate figures

        if figno == 1
            model.options.IMPSELECT = [1,2,18,21,29,30,40,15,33]; %[15,31,32,33,40,38,34,3]; 
            model.options.SHOCKSELECT = 1;
        else
            model.options.IMPSELECT = [1,2,18,21,29,30,40,15,33]; %[15,31,32,33,40,38,34,3];
            model.options.SHOCKSELECT = 3;
        end   

        if DO_IMPULSE1 == 1
            disp(' ')
            disp('STEP 7(a): Compute and plot Impulse responses')
            disp(' ')
            [ IRMAT ] = impulse_linear(model);
        end
        
        % SPECIAL PLOT FOR PAPER ... 
%         if sum(model.options.IMPSELECT - [15,40,33,31,18,29]) == 0
%             % Then do the RER decomposition IRF plot ...
%             
%             s_count = model.options.SHOCKSELECT;
%             s_count = ny+(nx-nshock)+s_count;
%             var_name = model.out.VARNAMES(model.options.IMPSELECT,:);
%             var_choose_current = [model.options.IMPSELECT,s_count];
%             createfigure_rer((0:model.options.HORIZON)', ...
%                                 IRMAT(:,var_choose_current),...
%                                     var_name, ...
%                                     model.out.VARNAMES(s_count,:),...
%                                     protocol);
%         end

    end
    
    
% % 7(b) Autocovariograms:
% % -----------------------------------------------------------------------

if MOMENTS == 1
        disp(' ')
        disp('STEP 7(b): Compute moments')
        disp(' ')
        
    varshock = blkdiag(zeros(nx-nshock,nx-nshock), OMEGA);
    
    % % Options:
    model.options.method = 0; % 0 == algebraic inversion; 
                              % 1 == doubling algorithm
    model.options.J = 1;
    
    model.options.MOM_SELECT_y = [15,33,38,3,36,34,40];  
    model.options.MOM_SELECT_x = []; 
    
    model.options.Y_INDEX = 36;
    model.options.NX1_INDEX = 26;
    model.options.NX2_INDEX = 27;
    
    % Cross-country correlation selector matrix: corr(x, x*)
    model.options.CROSSCOUNTRY = [  24, 26;         % (RER, NX)
                                    33, 25;         % (RER_CPI, C/C*)
                                    33, 15;         % (RER_CPI, E) 
                                    33, 26 ];       % (RER_CPI, NX)

    % Generate moments tables:
    [STD_table,XCF_Y_table,modelcorrcc,modelautocorrY] = ...
                               mom_out(model,DATA,varshock,GLOBAL_DISPLAY);

end

disp('---------PB_CALIBRA_RUN: THE END ----------------------------')

    