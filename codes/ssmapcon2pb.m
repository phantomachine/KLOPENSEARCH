function [lconval, nonlconval] = ssmapcon2pb(x,ALPHA,DELTA,SIGMA,...
                                             BETA,GAMA,THETA,KAPPA,...
                                             ZBAR,EPSILON,C,ETA,...
                                             OMEGA_I, OMEGA_F, ...
                                             TAU_X, TAU_H, TAU_K, b,...
                                             ~,~,~,~)
    
    % --------------------------------------------------------------------
    % SSMAPCON2PB.M
    % --------------------------------------------------------------------
    % Function solving for calibrated parameters in the proportional
    % bargaining model, given targets. 
    % Nonlinear equality constraints for use with fmincon.
    % --------------------------------------------------------------------
    %
    % INPUT:
    %   x : Row vector inputs
    %   PARAMETERS ...
    %
    % OUTPUT:
    %   lcoval     : Linear equality constraint vector value
    %   nonlconcal : Non-Linear equality constraint vector value
    %
    % --------------------------------------------------------------------
    % (c) 2010- T.Kam; Email: mortheus@gmail.com
    % 
    % See also SSMAPOBJ, F, F_H, F_K, G, C_Q, C_K, FMINCON
    
    % Extract candidate roots (minimizers)- check bounds in FMINCON usage:
    k = x(1);        % > 0
    
    A = x(2);        % > 0
    PHI = x(3);      % in (0,1)
    %SIGMA = x(4);    % in (0,1/2)
    B = x(4);        % > 0
    %ALPHA = x(6);    % in (0,1)
    THETA_B = x(5);  % in [0,1]
    H = x(6);        % in (0,1)
    KYratio = x(7);  % > 0
	%MDE = x(8);      % < 0
    NTS = x(8);     % in (0,1)
    MKP = x(9);     % in (0,m), m < +infty
    
    % (Ph*phi,K,q_czech,q,X) are functions of kss1:

    [~,K1,qz1,q1,X1,~,...
        ~,~,~,~,~,~,~,~,~,~,...
                        ~,~,~,~,~,markup1] ...
                                 = ssmapstatic_propbarg(OMEGA_I,OMEGA_F,...
                                                 EPSILON,THETA,...
                                                 TAU_X,TAU_H,TAU_K,...
                                                 A,B,C,...
                                                 ALPHA,ZBAR,DELTA,PHI,...
                                                 ETA,GAMA,SIGMA,KAPPA,...
                                                 BETA,THETA_B,b,k);
    
    % Real GDP definition in terms of K/Y ratio ... k
    
    
    REL_PRICE = (OMEGA_I^((EPSILON-1)/EPSILON))*THETA; % Ph/P
    
    s_Kinv = REL_PRICE*ZBAR*(K1/H)^(ALPHA-1)...
                  *( 1 + ((SIGMA*(1-TAU_H))/A)...
                  *K1*(KAPPA*c_q(q1,K1,PHI)+(1-KAPPA)*c_q(qz1,K1,PHI))/K1 );       
                    
    f1 = 1/KYratio - s_Kinv;
   
    % Euler equation ... PHI
    f2 = -DELTA + (1-1/BETA)/(1-TAU_K) ...
          + REL_PRICE*ZBAR*F_k(K1,H,ALPHA)...
          - (SIGMA*(1+TAU_X)/((1-TAU_K)*U_X(X1,B,GAMA))) ...
            *(KAPPA*c_k(q1,K1,PHI) + (1-KAPPA)*c_k(qz1,K1,PHI));
        
    % Resource constraint and y_h market clearing ... A
    f3 = -A + (1+TAU_X)^(-1) *(1-TAU_H)*B*REL_PRICE*ZBAR*F_h(K1,H,ALPHA)...
            *(((1+TAU_X)/H)^GAMA)...
         *((OMEGA_I/OMEGA_F - ((1-ALPHA)*TAU_H+ALPHA*TAU_K)*REL_PRICE)...
            *ZBAR*F(K1,H,ALPHA)/H - (1-TAU_K)*DELTA*(K1/H))^(-GAMA);
    
    
      
   % DM "nontradables" consumption share in total DM+CM consumption ... B
      
        qDM = ((SIGMA*(1-TAU_H))/A)...
             *REL_PRICE*ZBAR*((K1/H)^ALPHA)...
              *(KAPPA*c_q(q1,K1,PHI)*q1+(1-KAPPA)*c_q(qz1,K1,PHI)*qz1);
              
    f5 = NTS - qDM/(X1 + qDM);
    
   
    
    % Markup
    
    f6 = MKP - markup1;
    
    % Collect equations
    lconval = [];                          % Linear binding constraints
    nonlconval = [f1,f2,f3,f5,f6];   % Nonlinear binding constraints