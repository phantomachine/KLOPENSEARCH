% PT_SS

settings;

DisplayFMINCON = 'iter-detailed'; % {'iter', 'final','off', 'notify',...
%     'iter-detailed', 'final-detailed'}

% % COMPECON_SS.m
% % -----------------------------------------------------------------------
% % This script performs calibration and computes the RCE's steady state
% % allocations.
% %
% % -----------------------------------------------------------------------
% % (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------
% %
% % See also COMPECON_RUN, COMPECON_EQ, PARETO_RUN, PARETO_EQ, PARETO_SS,
% %          SSMAPOBJ2, SSMAPCON2, SSMAPSTATIC, FMINCON

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('COMPECON_SS.M: STEP 1. SET MEASURABLE PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')

% Invoke settings in PARAMETERSCOMMON.M:
parameterscommon;

% Special case (Table 4, Column (ii) NO ANONYMITY)
%KAPPA = 0.0001;

% Estimate ALPHA and SIGMA:
ALPHA = 0.33;   % Fix to match labor share of 2/3
SIGMA = 0.078;

H_target = 0.28;

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('COMPECON_SS.M: STEP 2. CALIBRATING NON-MEASURABLE PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')


options0 = optimset('Display', DisplayFMINCON,...
    'Algorithm','active-set',...
    'TolFun', 1e-15, ...
    'TolCon',1e-15, ...
    'TolX',1e-15, ...
    'MaxFunEvals',5e+4, ...
    'MaxIter',5e+4);

% Name calibration parameters and targets:

NAME_cabparams = {'k','A','PHI','B'};
N_par = length(NAME_cabparams);

NAME_cabtarget = {'H','K/Y','NTS'};
N_tar = length(NAME_cabtarget);

% Feasible space for minimizers:
% {'k','A','PHI','SIGMA','ALPHA','B'}, {'H','K/Y','v','LS','NTS'}

% Feasible parameter space: Lower bounds
k_lb        = 1e-12;
A_lb        = 1e-12;
PHI_lb      = 1.05;%1.0005; %1.05;
%         SIGMA_lb    = 1e-12;
B_lb        = 1e-12;
H_lb        = 1e-12;
KY_lb       = 1e-12;
NTS_lb      = 1e-12;

lb = [k_lb,A_lb,PHI_lb,B_lb,...
    H_lb,KY_lb,NTS_lb];

% Feasible parameter space: Upper bounds
k_ub        = 100;
A_ub        = 100;
PHI_ub      = 1.6; %1.001; %1.6;
%         SIGMA_ub    = 0.26/4;
B_ub        = 10;
H_ub        = 0.4;
KY_ub       = 15;
NTS_ub      = 0.999;

ub = [k_ub,A_ub,PHI_ub,B_ub,...
                        H_ub,KY_ub,NTS_ub];

% initial guesses

x0 = [10,14,1.1,1,0.3,2,0.1];

% Do constraint minimization - calibration (1st moment matching)
[xss,fmin,exitflag,fmcoutput,lm] = fmincon(@ssmapobj1pt,x0,[],[],[],[],...
                                            lb,ub,@ssmapcon1pt,options0,...
                                            ALPHA,SIGMA,DELTA,BETA,...
                                            GAMA,THETA,KAPPA,...
                                            ZBAR,EPSILON,C,ETA,...
                                            OMEGA_I, OMEGA_F, ...
                                            TAU_X, TAU_H, TAU_K, ...
                                            H_target,KY_target, ...
                                            NTS_target);

cabparams = xss(1 : N_par);
cabtarvar = xss(end-N_tar+1 : end);
cabtarget = [H_target, KY_target,NTS_target];

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
%SIGMA = xss(4);
%ALPHA = xss(5);
B = xss(4);

H1 = xss(5);
KY1 = xss(6);

% (Ph*phi,K,q_czech,q,X) are functions of kss1:

[REL_Pyh1,K1,qz1,q1,X1,Yh1,...
    I1,yh1,yf1,yf2,yh2,REL_Pyf1,w1,RER,E,P1,...
    Pyh1,Pyf1,Pyf2,Pyh2,NX1] ...
                                   = ssmapstaticpt(OMEGA_I,OMEGA_F,...
                                                    EPSILON,THETA,...
                                                    TAU_X,TAU_H,TAU_K,...
                                                    A,B,C,...
                                                    ALPHA,ZBAR,DELTA,PHI,...
                                                    ETA,GAMA,SIGMA,KAPPA,...
                                                    BETA,k1,H1);


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

RELX = X1/X2;

TOT = Pyf1/Pyh1;

Ptilde_q1 = ((A*psi1/(P1*(1-TAU_H)*w1))^(-1))*c_q(q1,K1,PHI);
% Home money trades: price
% NOTE: c_q, Marginal Cost in utils
Ptilde_qz1 = ((A*psi1/(P1*(1-TAU_H)*w1))^(-1))*c_q(qz1,K1,PHI);

Ptilde1 = ( KAPPA*Ptilde_q1 + (1-KAPPA)*Ptilde_qz1 );
% Average DM price level: Home

Ptilde_q2 = ((A*psi2/(P2*(1-TAU_H)*w2))^(-1))*c_q(q2,K2,PHI);
% Foreign money trades: price
Ptilde_qz2 = ((A*psi2/(P2*(1-TAU_H)*w2))^(-1))*c_q(qz2,K2,PHI);
% Foreign credit trades: price

Ptilde2 = ( KAPPA*Ptilde_q2 + (1-KAPPA)*Ptilde_qz2 );

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

% Average DM price level: Foreign

RER_YPI = E* YPI2 / YPI1;

Agg_C1 = ((P1*X1+SIGMA*(1-KAPPA)*Ptilde_qz1*qz1+...
                          SIGMA*KAPPA*Ptilde_q1*q1)/YPI1);
Agg_C2 = ((P2*X2+SIGMA*(1-KAPPA)*Ptilde_qz2*qz2+...
                          SIGMA*KAPPA*Ptilde_q2*q2)/YPI2);
Agg_Y1 = (Pyh1*yh1+SIGMA*(KAPPA*Ptilde_q1*q1+(1-KAPPA)*Ptilde_qz1*qz1))/YPI1;
Agg_Y2 = (Pyf2*yf2+SIGMA*(KAPPA*Ptilde_q2*q2+(1-KAPPA)*Ptilde_qz2*qz2))/YPI2;

Agg_I1 = P1*I1/YPI1;
Agg_I2 = P2*I2/YPI2;

P2P1 = YPI2/YPI1;

% DM Hours/effort
    H_dm1 = SIGMA*(KAPPA*cost(q1,K1,PHI) + (1-KAPPA)*cost(qz1,K1,PHI));
    
    % Total Labor hours
    HTOT1 = H1 + H_dm1;
    HTOT2 = HTOT1; 
    
% MDE: Money demand interest elasticity
THETA_B = 0;

[m_elasticity] = m_interest_elasticity(q1,qz1,K1,H1,ZBAR,X1,...
                                            ALPHA,OMEGA_I,OMEGA_F,...
                                            THETA,DELTA,...
                                            EPSILON,THETA_B,...
                                            ETA,A,b,B,C,SIGMA,KAPPA,PHI,...
                                            GAMA,...
                                            TAU_X,TAU_H,TAU_K,...
                                            'pt');
                      

%disp('Total Labor income share')
H_Income_Total = P1*ALPHA*REL_Pyh1*FF1 ...    
                            + SIGMA*(Ptilde_q1* KAPPA* q1*(1/PHI)) ...
                            + Ptilde_qz1*(1-KAPPA)*qz1*(1-1/PHI);

Labor_Share = H_Income_Total/(Agg_Y1*YPI1);


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
            m_elasticity, m_elasticity;
            Labor_Share, Labor_Share    ];

N_SSMAT = size(SSMAT,1);

NAME_SSMAT = {'K','X','H','HTOT','q_cz','q','I','w',...
    'P','Py_ii','Py_ji','y_ii','y_ji','NX','GDP','MDE','LS'};

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

% ACTIVATE: For Exercise in Table 4(iii)- Anonimity + No DM complementarity
%PHI = 1.0000001;

