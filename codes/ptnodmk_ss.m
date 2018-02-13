% PTNODMK_SS.M

settings;

DisplayFMINCON = 'iter-detailed'; % {'iter', 'final','off', 'notify',...
                                   %     'iter-detailed', 'final-detailed'}

% % PTNODMK_SS.m
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
disp('PTNODMK_SS.M: STEP 1. SET BASELINE PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')

    % Invoke settings in PARAMETERSCOMMON.M:
    parameterscommon;
  
    
    % Calibration from pt_ss ...
    % (Pricetaking baseline with credit, money and DM capital):
        ALPHA = 0.33;
        A = 0.3326;
        B = 0.0967;
        SIGMA = 0.078;
        %PHI = 1.289;   
        PHI = 1.000000000001;


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('PTNODMK_SS.M: STEP 2. SOLVE FOR STEADY STATE ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')


% Solve Euler Equation for k := K/H
options0 = optimset( 'Display','off',...
                     'TolFun',10e-12,...
                     'MaxFunEvals',10e+5,...
                     'MaxIter',10e+5         );

k_init = 2;
                
[kss, fss, flag, model.fsout] = fsolve('keulerpt',k_init,options0,...
                                    A,B,C,ETA,SIGMA,KAPPA,...
                                    ALPHA,BETA,GAMA,DELTA,ZBAR,EPSILON,...
                                    THETA,OMEGA_I,OMEGA_F,PHI,...
                                    TAU_X,TAU_K,TAU_H);

if flag && strcmp(options0.Display,'off')   
    fprintf('FSOLVE: Equation Solved.\n\t Type model.fsout for more...\n');  
else 
    warning('Error:FSOLVE','You may not have true zero')
end

    
    % Pin down other steady state variables:

    k1 = kss;

    % (Ph*phi,K,q_czech,q,X) are functions of kss1:

    [REL_Pyh1,K1,qz1,q1,X1,Yh1,...
        I1,yh1,yf1,yf2,yh2,REL_Pyf1,w1,RER,E,P1,...
                        Pyh1,Pyf1,Pyf2,Pyh2,NX1] ...
                                   = ssmapstaticptnodmk(OMEGA_I,OMEGA_F,...
                                                 EPSILON,THETA,...
                                                 TAU_X,TAU_H,TAU_K,...
                                                 A,B,C,...
                                                 ALPHA,ZBAR,DELTA,PHI,...
                                                 ETA,GAMA,SIGMA,KAPPA,...
                                                 BETA,k1);

H1 = k1*K1;

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

% DM Hours/effort
    H_dm1 = SIGMA*(KAPPA*cost(q1,K1,PHI) + (1-KAPPA)*cost(qz1,K1,PHI));
    
    % Total Labor hours
    HTOT1 = H1 + H_dm1;
    HTOT2 = HTOT1; 
                      

% Check for symmetry
SSMAT = [   K1,       K2;
            X1,       X2;
            H1,       H2;
            HTOT1,    HTOT2;
            qz1,      qz2;
            q1,       q2;
            I1,       I2;
            w1,       w2;
            P1,       P2;
            Pyh1,   Pyf2;
            Pyf1,   Pyh2;
            yh1,     yf2;
            yf1,     yh2;
            NX1,     NX2;
            Agg_Y1, Agg_Y2    ];

N_SSMAT = size(SSMAT,1);

NAME_SSMAT = {'K','X','H','HTOT','q_cz','q','I','w',...
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

THETA_B = 0;

[m_elasticity] = m_interest_elasticity(q1,qz1,K1,H1,ZBAR,X1,...
                                            ALPHA,OMEGA_I,OMEGA_F,...
                                            THETA,DELTA,...
                                            EPSILON,THETA_B,...
                                            ETA,A,b,B,C,SIGMA,KAPPA,PHI,...
                                            GAMA,...
                                            TAU_X,TAU_H,TAU_K,...
                                            'pt');
                                        
disp('Money interest semi-elasticity = ')
disp(m_elasticity)

