function [REL_Pyh1,K1,q_czech1,q1,X1,Yh1,I1,yh1,yf1,yf2,yh2,...
            REL_Pyf1,w1,RER,E,P1,Pyh1,Pyf1,Pyf2,Pyh2,NX1] = ...
                                     ssmapstaticptnodmk(OMEGA_I,OMEGA_F,...
                                                 EPSILON,THETA,...
                                                 TAU_X,TAU_H,TAU_K,...
                                                 A,B,C,...
                                                 ALPHA,ZBAR,DELTA,PHI,...
                                                 ETA,GAMA,SIGMA,KAPPA,...
                                                 BETA,k)

    % Composite parameters:
    % ---------------------------------------------------------------------
    
    REL_Pyh = (OMEGA_I^((EPSILON-1)/EPSILON))*THETA; % Ph/P
    
    % Side equations to pin down (q, q_czech, K, X):
    
    K = (((1+TAU_X)*(1-TAU_H)*B*REL_Pyh/A)*ZBAR*k^ALPHA)^(1/GAMA)...
          /(((1+TAU_X)^(-1))...
            *(OMEGA_I/OMEGA_F-((1-ALPHA)*TAU_H+ALPHA*TAU_K)*REL_Pyh)...
            *ZBAR*k^(ALPHA-1) - (1-TAU_K)*DELTA);
    
    q_czech = ((C*K^(PHI-1))/PHI)^(1/(ETA+PHI-1));
    
    q = ((SIGMA*KAPPA*C*K^(PHI-1))...
                    /(PHI*(1/BETA-1+SIGMA*KAPPA)) )^(1/(ETA+PHI-1));
                
    X = ((1+TAU_X)*(1-TAU_H)*B*REL_Pyh*ZBAR*k^ALPHA/A)^(1/GAMA);
    
    H = K/k;
    
    % Other steady state equations: Symmetric global SME
    % ---------------------------------------------------------------------
    
    % Just relabelling:
    REL_Pyh1 = REL_Pyh;
    K1 = K;
    q_czech1 = q_czech;
    q1 = q;
    X1 = X;
    
    
    Yh1 = ZBAR*F(K,H,ALPHA); % CM output
     
    I1 = DELTA*K; % Investment
    
    yh1 = Yh1/OMEGA_F; % Home-produced indetermediate good demand
    
    yf1 = (THETA/(1-THETA))^(EPSILON/(1-EPSILON))*yh1;
                     % Foreign-produced indetermediate good demand
                      
    REL_Pyf1 = GG_2(yh1,yf1,EPSILON,THETA);    % Pf1/P       
    
    w1 = REL_Pyh1*ZBAR*F_h(K,H,ALPHA);        % MP-labor
    
    epsi1 = 1; % Money supply growth
    
    P1 = A*epsi1 /(w1*(1-TAU_H)*g(q,K,ZBAR,PHI));      % inverse real balance
    P2 = P1;
    
    RER = 1;       % RER = U_X(X*)/U_X(X) = 1, since X* = X
    E = 1;
    
    Pyh1 = REL_Pyh1*P1;  % P_yh, definitional
    Pyf1 = REL_Pyf1*P1;  % P_yf, definitional         
    
    yh2 = yf1; % Symmetry in intermediate goods demand
    yf2 = yh1;
    
    REL_Pyf2 = GG_1(yf2,yh2,EPSILON,THETA);      % P_yf2*/P*
    REL_Pyh2 = GG_2(yf2,yh2,EPSILON,THETA);      % P_yh2*/P*
    
    Pyf2 = REL_Pyf2*P2;
    Pyh2 = REL_Pyh2*P2;
    
    NX1 = Pyh2*E*yh2/P1 - Pyf1*yf1/P1; % i.e. NX = 0
    
    %GDP1 = REL_Pyf1*Yh1 + NX1; % GDP
    
    % % DM markup
    % % markup = 0 for price taking
    % % markup > 0 for bargaining
    % 
    % markup = g(q1,K1,PHI)/(q1*c_q(q1,K1,PHI)) - 1; 