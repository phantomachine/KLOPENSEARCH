settings;

DisplayFMINCON = 'iter-detailed'; % {'iter', 'final','off', 'notify',...
                                   %     'iter-detailed', 'final-detailed'}

% % COMPECONCREDITPROPB_SS.m
% % -----------------------------------------------------------------------
% % This script performs calibration and computes the RCE's steady state 
% % allocations. Model with proportional bargaining
% % 
% % -----------------------------------------------------------------------
% % (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------
% %
% % See also COMPECONPROPB_RUN, COMPECONPROPB_EQ, SSMAPOBJ3PROPBARG, 
% %          SSMAPCON3PROPBARG, SSMAPSTATICPROPBARG, FMINCON 

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('COMPECON_SS.M: STEP 1. SET MEASURABLE PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')

    % Invoke settings in PARAMETERSCOMMON.M:
    parameterscommon;

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('COMPECON_SS.M: STEP 2. CALIBRATING NON-MEASURABLE PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')


    options0 = optimset('Display', DisplayFMINCON,...
					'Algorithm','active-set',...
                                        'TolFun', 1e-8, ...
                                        'TolCon',1e-8, ...
                                        'TolX',1e-8, ...
                                        'MaxFunEvals',5e+4, ...
                                        'MaxIter',5e+10); 

    % Name calibration parameters and targets:

    NAME_cabparams = {'k','A','PHI','SIGMA','l','B','BTHETA'}; 
    N_par = length(NAME_cabparams);

    NAME_cabtarget = {'H','K/Y','v','LS','NTS','MUP'};      
    N_tar = length(NAME_cabtarget);

    % Feasible space for minimizers:
    % {'k','A','PHI','SIGMA','ALPHA','B'}, {'H','K/Y','v','LS','NTS'}

        % Feasible parameter space: Lower bounds
        k_lb        = 1e-12;
        A_lb        = 1e-12;
        PHI_lb      = 1.05;
        SIGMA_lb    = 1e-12;
        ALPHA_lb    = 1e-12;
        B_lb        = 1e-12;
        BTHETA_lb   = 1e-12;
        H_lb        = 1e-12;
        KY_lb       = 1e-12;
        v_lb        = 1e-12;
        LS_lb       = 1e-12;
        NTS_lb      = 1e-12;
        MUP_lb      = 1e-12;

        lb = [k_lb,A_lb,PHI_lb,SIGMA_lb,ALPHA_lb,B_lb,BTHETA_lb,...
                                       H_lb,KY_lb,v_lb,LS_lb,NTS_lb,MUP_lb];

        % Feasible parameter space: Upper bounds
        k_ub        = 10;%100;
        A_ub        = 100;
        PHI_ub      = 5;
        SIGMA_ub    = 0.4999;
        ALPHA_ub    = 0.999;
        B_ub        = 10;
        BTHETA_ub   = 0.999;
        H_ub        = 0.4;
        KY_ub       = 6;
        v_ub        = 10;
        LS_ub       = 0.8;
        NTS_ub      = 0.999;
        MUP_ub      = 5;
        
        ub = [k_ub,A_ub,PHI_ub,SIGMA_ub,ALPHA_ub,B_ub,BTHETA_ub,...
                                    H_ub,KY_ub,v_ub,LS_ub,NTS_ub,MUP_ub];

    % initial guesses
    k0 = 5;%10;
    A0 = 14;
    PHI0 = 1.5;%1.1;
    SIG0 = 0.26/4;
    ALPHA0 = 0.288;
    B0 = 1;
    BTHETA0 = 0.95;
    H0 = 0.3;
    KY0 = 2;
    v0 = 1;
    LS0 = 0.5;
    NTS0 = 0.4;
    MUP0 = 0.3;
    
    x0 = [k0, A0, PHI0, SIG0, ALPHA0, B0, BTHETA0, ...
                                        H0, KY0, v0, LS0, NTS0, MUP0];
    
    % Do constraint minimization - calibration (1st moment matching)
    [xss,fmin,exitflag,fmcoutput,lm] = fmincon(@ssmapobj3propbarg,x0,[],[],[],[],...
                                            lb,ub,@ssmapcon3propbarg,options0,...
                                            DELTA,BETA,GAMA,THETA,KAPPA,...
                                                ZBAR,EPSILON,C,ETA,...
                                                b, ...
                                                OMEGA_I, OMEGA_F, ...
                                                TAU_X, TAU_H, TAU_K, ...
                                                   H_target, KY_target, ...
                                            v_target,LS_target,...
                                            NTS_target,MUP_target);

    cabparams = xss(1 : N_par);
    cabtarvar = xss(end-N_tar+1 : end);
    cabtarget = [H_target, KY_target,v_target,LS_target,NTS_target,MUP_target];

    if exitflag ~= 0 || exitflag ~= -1 || exitflag ~= -2
        if strcmp(GLOBAL_DISPLAY, 'on') == 1
            disp(' ')

            disp('        Table 1: Optimized Calibrated Parameters')
            disp('-------------------------------------------------------')

            for i = 1:N_par
                disp([NAME_cabparams{i}, sprintf('\t\t %-6.3f',...
                                                            cabparams(i))])
            end

            disp('-------------------------------------------------------')
            disp('Note: k is endogenously determined.')
            disp('=======================================================')
            disp(' ')

            disp('        Table 2: Calibration Target Variables')
            disp('-------------------------------------------------------')
            disp(['Variables  ', 'Optimized Values  ', 'Empirical Target'])
            disp('-------------------------------------------------------')
            for i = 1:N_tar
                disp([NAME_cabtarget{i},...
                            sprintf('\t\t %-6.3f',cabtarvar(i)),...
                                      sprintf('\t\t %-6.3f',cabtarget(i))])
            end
            disp('=======================================================')
        end
    elseif exitflag == 0 || exitflag == -1 || exitflag == -2
        warning('COMPECON_SS:Error','Calibration task not completed...')
        break
    end

    % Pin down other steady state variables:

    k1 = xss(1);
    A = xss(2);
    PHI = xss(3);
    SIGMA = xss(4);
    ALPHA = xss(5);
    B = xss(6);
    BTHETA = xss(7);

    H1 = xss(8);
    KY1 = xss(9);

    % (Ph*phi,K,q_czech,q,X) are functions of kss1:

    [REL_Pyh1,K1,qz1,q1,X1,Yh1,...
        I1,yh1,yf1,yf2,yh2,REL_Pyf1,w1,RER,E,P1,...
                        Pyh1,Pyf1,Pyf2,Pyh2,NX1,GDP1] ...
                                   = ssmapstaticpropbarg(OMEGA_I,OMEGA_F,...
                                                 EPSILON,THETA,...
                                                 TAU_X,TAU_H,TAU_K,...
                                                 A,B,C,...
                                                 ALPHA,ZBAR,DELTA,PHI,...
                                                 ETA,GAMA,SIGMA,KAPPA,...
                                                 BETA,BTHETA,k1,H1,KY1);

   
    % By symmetry, for Foreign's steady state SME allocation:
    k2 = k1;
    X2 = X1;
    K2 = K1;
    H2 = H1;
    qz2 = qz1;
    q2 = q1;
    
    Z1 = ZBAR;
    Z2 = Z1;
    
    psi1 = EPSIBAR;
    psi2 = 1;
    
    I2 = I1;
    w2 = w1;
    P2 = P1;
    NX2 = NX1;
    
    Y1 = GDP1;
    Y2 = Y1;
    
    RELX = X1/X2;
    
    TOT = Pyf1/Pyh1;
    
    % Check for symmetry
    SSMAT = [ K1,       K2;
              X1,       X2;
              H1,       H2;
              qz1,     qz2;
              q1,       q2;
              I1,       I2;
              w1,       w2;
              P1,       P2;
              Pyh1,   Pyf2;
              Pyf1,   Pyh2;
              yh1,     yf2;
              yf1,     yh2;
              NX1,     NX2;
              Y1,      Y2 ];
          
              N_SSMAT = size(SSMAT,1);
    
    NAME_SSMAT = {'K','X','H','q_cz','q','I','w',...
                         'P','Py_ii','Py_ji','y_ii','y_ji','NX','GDP'};
    
    if strcmp(GLOBAL_DISPLAY, 'on') == 1
        disp(' ')
        disp('            Table 3: Steady State Values')
        disp('-------------------------------------------------------')
        disp(['Variables ', ' Home (i=h,j=f)   ', 'Foreign (i=f,j=h)'])
        disp('-------------------------------------------------------')
        for i = 1:N_SSMAT
            disp([NAME_SSMAT{i},sprintf('\t\t %-6.3f',SSMAT(i,1)),...
                                    sprintf('\t\t %-6.3f',SSMAT(i,2))])
        end
        disp('=======================================================')
    end