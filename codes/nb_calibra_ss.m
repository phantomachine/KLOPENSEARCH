
settings;

DisplayFMINCON = 'iter-detailed'; % {'iter', 'final','off', 'notify',...
                                   %     'iter-detailed', 'final-detailed'}

% % % NB_CALIBRA_SS.M
% % -----------------------------------------------------------------------
% % This script performs calibration and computes the RCE's steady state 
% % allocations.
% % 
% % -----------------------------------------------------------------------
% % (c) 2010 - , Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------
% %
% % See also NB_CALIBRA_RUN, NB_EQ, SSMAPSTATIC_NASHBARG, FMINCON 

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('NB_CALIBRA_SS.M: STEP 1. SET BASELINE PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')

    % Invoke settings in PARAMETERSCOMMON.M:
    parameterscommon;
  
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('NB_CALIBRA_SS.M: STEP 2. SOLVE FOR STEADY STATE ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')


options0 = optimset('Display', DisplayFMINCON,...
                    'TolFun', 1e-8, ...
                    'TolCon',1e-8, ...
                    'TolX',1e-8, ...
                    'Algorithm','sqp', ...
                    'MaxFunEvals',5e+4, ...
                    'MaxIter',5e+4);

% The magic hand ...                
ALPHA = 0.3;
SIGMA = 0.03;

H_target = 0.3;

% Name calibration parameters and targets:

NAME_cabparams = {'k','q','A','PHI','B','THETA_B'};
N_par = length(NAME_cabparams);

NAME_cabtarget = {'H','K/Y','NTS','MKP'};
N_tar = length(NAME_cabtarget);

% Feasible parameter space: Lower bounds
k_lb        = 1e-12;
q_lb        = 0.001;
A_lb        = 1e-12;
PHI_lb      = 1.05;
B_lb        = 1e-12;
THETA_B_lb  = 0.95;
H_lb        = 1e-12;
KY_lb       = 1e-12;
NTS_lb      = 0.45;
MKP_lb      = 0.3;


lb = [k_lb,q_lb,A_lb,PHI_lb,B_lb,THETA_B_lb,...
                                   H_lb,KY_lb,NTS_lb,MKP_lb];

% Feasible parameter space: Upper bounds
k_ub        = 100;
q_ub        = 1.5;
A_ub        = 100;
PHI_ub      = 1.6;
B_ub        = 10;
THETA_B_ub  = 1;
H_ub        = 0.3;
KY_ub       = 8.92;
NTS_ub      = 0.52;
MKP_ub      = 0.35;

ub = [k_ub,q_ub,A_ub,PHI_ub,B_ub,THETA_B_lb,...
                                  H_ub,KY_ub,NTS_ub,MKP_ub];

% initial guesses
k0 = 16.684;
q0 = 0.263;
A0 = 0.372;
PHI0 = 1.277;
THETA_B0 = 0.8999;
B0 = 0.119;
H0 = 0.3;
KY0 = 8.8;
NTS0 = 0.48;
MKP0 = 0.3;

x0 = [k0, q0, A0, PHI0, THETA_B0, B0, H0, KY0, NTS0, MKP0 ];

% Do constraint minimization - calibration (1st moment matching)
[xss,fmin,exitflag,fmcoutput,lm] = fmincon(@ssmapobj2nb,x0,...
                                            [],[],[],[],...
                                            lb,ub,...
                                            @ssmapcon2nb,...
                                            options0,...
                                            ALPHA,DELTA,SIGMA,BETA,...
                                            GAMA,THETA,KAPPA,...
                                            ZBAR,EPSILON,C,ETA,...
                                            OMEGA_I, OMEGA_F, ...
                                            TAU_X, TAU_H, TAU_K, b, ...
                                            H_target,KY_target, ...
                                            NTS_target,...
                                            MKP_target);

cabparams = xss(1 : N_par);
cabtarvar = xss(end-N_tar+1 : end);
cabtarget = [H_target, KY_target,NTS_target,MKP_target];

if exitflag ~= 0 || exitflag ~= -1 || exitflag ~= -2
    if strcmp(GLOBAL_DISPLAY, 'on') == 1
        disp(' ')
        
        disp('        Table 1: Optimized Calibrated Parameters       ')
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
    warning('NB_CALIBRA_SS:Error','Calibration task not completed...')
    break
end

% Pin down other steady state variables:

    k1 = xss(1);        % > 0
    q1 = xss(2);        % > 0
    A = xss(3);        % > 0
    PHI = xss(4);      % in (0,1)
    B = xss(5);        % > 0
    THETA_B = xss(6);  % in [0,1]
    H1 = xss(7);        % in (0,1)
    KYratio = xss(8);  % > 0
	NTS = xss(9);     % in (0,1)
    MKP = xss(10);     % in (0,m), m < +infty

[REL_Pyh1,K1,qz1,X1,Yh1,I1,yh1,yf1,yf2,yh2,...
            REL_Pyf1,w1,RER,E,P1,Pyh1,Pyf1,Pyf2,Pyh2,NX1,~] ...
                                 = ssmapstatic_nashbarg(OMEGA_I,OMEGA_F,...
                                                 EPSILON,THETA,...
                                                 TAU_X,TAU_H,TAU_K,...
                                                 A,B,C,...
                                                 ALPHA,ZBAR,DELTA,PHI,...
                                                 ETA,GAMA,SIGMA,KAPPA,...
                                                 BETA,THETA_B,b,k1,q1);
    
    % DM Hours/effort
    H_dm1 = SIGMA*(KAPPA*cost(q1,K1,PHI) + (1-KAPPA)*cost(qz1,K1,PHI));
    
    % Total Labor hours
    HTOT1 = H1 + H_dm1;
    HTOT2 = HTOT1; 
    
    markup1 = MKP;
    markup2 = markup1;
    
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
    
    %Y1 = GDP1;
    %Y2 = Y1;
    
    RELX = X1/X2;
    
    TOT = Pyf1/Pyh1;
   
    Ptilde_q1 = ((A*psi1/(P1*(1-TAU_H)*w1))^(-1))...
                            *g_prop(q1,K1,Z1,THETA_B,PHI,ETA,b,C)/q1;     
                                 % Home money trades: price
                                 % NOTE: c_q, Marginal Cost in utils
    Ptilde_qz1 = ((A*psi1/(P1*(1-TAU_H)*w1))^(-1))*c_q(qz1/Z1,K1,PHI)/Z1;    
    
    Ptilde1 = ( KAPPA*Ptilde_q1 + (1-KAPPA)*Ptilde_qz1 );
                                             % Average DM price level: Home
     
    Ptilde_q2 = ((A*psi2/(P2*(1-TAU_H)*w2))^(-1)) ...
                            *g_prop(q2,K2,Z2,THETA_B,PHI,ETA,b,C)/q2;
                                         % Foreign money trades: price
    Ptilde_qz2 = ((A*psi2/(P2*(1-TAU_H)*w2))^(-1))*c_q(qz2/Z2,K2,PHI)/Z2;   
                                         % Foreign credit trades: price
    
    Ptilde2 = ( KAPPA*Ptilde_q2 + (1-KAPPA)*Ptilde_qz2 );
                                         % Average DM price level: Foreign   
        
NTS_s1 = 1-Pyh1*yh1 ...
           /(Pyh1*yh1+SIGMA*(KAPPA*Ptilde_q1*q1+(1-KAPPA)*Ptilde_qz1*qz1));
NTS_s2 = 1-Pyf2*yf2 ...
           /(Pyf2*yf2+SIGMA*(KAPPA*Ptilde_q2*q2+(1-KAPPA)*Ptilde_qz2*qz2));

% Intermediate goods production output
FF1 = Z1*F(K1,H1,ALPHA);
FF2 = Z2*F(K2,H2,ALPHA);

% Total output in CM units of consumption
Y1 = (Pyh1*FF1 - SIGMA*KAPPA*Ptilde_q1*q1...
                                    - SIGMA*(1-KAPPA)*Ptilde_qz1*qz1)/P1;
Y2 = (Pyf2*FF2 - SIGMA*KAPPA*Ptilde_q2*q2...
                                    - SIGMA*(1-KAPPA)*Ptilde_qz2*qz2)/P2;                                
                               

% GDP deflator: price levels
YPI1 = (1-NTS_s1)*P1 + NTS_s1*Ptilde1;
YPI2 = (1-NTS_s2)*P2 + NTS_s2*Ptilde2;

P2P1 = YPI2/YPI1;

% Average DM price level: Foreign

RER_YPI = E* YPI2 / YPI1;

Agg_C1 = ((P1*X1+SIGMA*(1-KAPPA)*Ptilde_qz1*qz1+...
                          SIGMA*KAPPA*Ptilde_q1*q1)/YPI1);
Agg_C2 = ((P2*X2+SIGMA*(1-KAPPA)*Ptilde_qz2*qz2+...
                          SIGMA*KAPPA*Ptilde_q2*q2)/YPI2);
Agg_Y1 = (Pyh1*yh1+SIGMA*(KAPPA*Ptilde_q1*q1+(1-KAPPA)*Ptilde_qz1*qz1))...
                                                                     /YPI1;
Agg_Y2 = (Pyf2*yf2+SIGMA*(KAPPA*Ptilde_q2*q2+(1-KAPPA)*Ptilde_qz2*qz2))...
                                                                     /YPI2;

Agg_I1 = P1*I1/YPI1;
Agg_I2 = P2*I2/YPI2;


% Total Dm and CM labor income share ... ALPHA
        
    H_Income_Total = P1*ALPHA*REL_Pyh1*FF1 ...    
                            + SIGMA*(Ptilde_q1* KAPPA* q1*(1/PHI)) ...
                            + Ptilde_qz1*(1-KAPPA)*qz1*(1-1/PHI);
    LS = H_Income_Total/(Agg_Y1*YPI1);
    

    
    % Real money demand interest elasticity ... SIGMA
          
    MDE = m_interest_elasticity(q1,qz1,K1,H1,ZBAR,X1,...
                                            ALPHA,OMEGA_I,OMEGA_F,...
                                            THETA,DELTA,...
                                            EPSILON,THETA_B,...
                                            ETA,A,b,B,C,SIGMA,KAPPA,PHI,...
                                            GAMA,...
                                            TAU_X,TAU_H,TAU_K,...
                                            'nb');
                   

% Check for symmetry
SSMAT = [   K1,       K2;
            X1,       X2;
            H1,       H2;
            HTOT1,    HTOT2;
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
            Agg_Y1, Agg_Y2;
            MKP,     MKP;
            MDE,     MDE;
            LS,      LS      ];

N_SSMAT = size(SSMAT,1);

NAME_SSMAT = {'K','X','H','HTOT','q_cz','q','I','w',...
    'P','Py_ii','Py_ji','y_ii','y_ji','NX','GDP', 'Markup','MDE','LS'};

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

[m_elasticity] = m_interest_elasticity(q1,qz1,K1,H1,ZBAR,X1,...
                                            ALPHA,OMEGA_I,OMEGA_F,...
                                            THETA,DELTA,...
                                            EPSILON,THETA_B,...
                                            ETA,A,b,B,C,SIGMA,KAPPA,PHI,...
                                            GAMA,...
                                            TAU_X,TAU_H,TAU_K,...
                                            'nb');
                                        