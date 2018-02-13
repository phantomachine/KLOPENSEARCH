function [lconval, nonlconval] = ssmapcon1pt(x,ALPHA,SIGMA,DELTA,BETA,GAMA,THETA,KAPPA,...
                                            ZBAR,EPSILON,C,ETA,...
                                            OMEGA_I, OMEGA_F, ...
                                            TAU_X, TAU_H, TAU_K,...
                                             ~,~,~)
    
    % --------------------------------------------------------------------
    % SSMAPCON.M
    % --------------------------------------------------------------------
    % Function solving for calibrated parameters in the price-taking
    % model, given targets. Nonlinear equality constraints
    % for use with fmincon.
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
    % See also SSMAPOBJ, F, F_H, F_K, G, C_Q, C_K
    
    % Extract candidate roots (minimizers)- check bounds in FMINCON usage:
    k = x(1);        % > 0
    
    A = x(2);        % > 0
    PHI = x(3);      % in (0,1)
%     SIGMA = x(4);    % in (0,1/2)
%     ALPHA = x(5);    % in (0,1)
    B = x(4);        % > 0
    
    H = x(5);        % in (0,1)
    KYratio = x(6);  % > 0
%     velocity = x(7); % > 0
%     LS = x(8);       % in (0,1)
    NTS = x(7);     % in (0,1)
    
    [REL_PRICE,K,q_czech,q,X,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] ...
                                   = ssmapstaticpt(OMEGA_I,OMEGA_F,...
                                                 EPSILON,THETA,...
                                                 TAU_X,TAU_H,TAU_K,...
                                                 A,B,C,...
                                                 ALPHA,ZBAR,DELTA,PHI,...
                                                 ETA,GAMA,SIGMA,KAPPA,...
                                                 BETA,k,H);
    
    % Real GDP definition in terms of K/Y ratio ... k
        s_Kinv = REL_PRICE*ZBAR*(K/H)^(ALPHA-1)...
                  *( 1 + ((SIGMA*(1-TAU_H))/A)...
                  *k*(KAPPA*c_q(q,K,PHI)+(1-KAPPA)*c_q(q_czech,K,PHI))/K );       
                    
    f1 = 1/KYratio - s_Kinv;
   
    % Euler equation ... PHI
    f2 = -DELTA + (1-1/BETA)/(1-TAU_K) ...
          + REL_PRICE*ZBAR*F_k(K,H,ALPHA)...
          - (SIGMA*(1+TAU_X)/((1-TAU_K)*U_X(X,B,GAMA))) ...
            *(KAPPA*c_k(q,K,PHI) + (1-KAPPA)*c_k(q_czech,K,PHI));
        
    % Resource constraint and y_h market clearing ... A
    f3 = -A + (1+TAU_X)^(-1) *(1-TAU_H)*B*REL_PRICE*ZBAR*F_h(K,H,ALPHA)...
            *(((1+TAU_X)/H)^GAMA)...
         *((OMEGA_I/OMEGA_F - ((1-ALPHA)*TAU_H+ALPHA*TAU_K)*REL_PRICE)...
            *ZBAR*F(K,H,ALPHA)/H - (1-TAU_K)*DELTA*(K/H))^(-GAMA);
    
   % Quantity theory of money equation ... SIGMA
          
          % Real GDP
   %       GDP = s_Kinv*K; 
      
          % Phi * M: from DM pricing FOC
   %       phi_M = ((1-TAU_H)*REL_PRICE*F_h(K,H,ALPHA)*g(q,K,PHI))/A;
      
   %   f4 = velocity - GDP / phi_M;
      
    % Labor share in CM output final goods terms: Cobb Douglas ... ALPHA
       
   %   f5 = ALPHA + log(LS*ZBAR)/log(K/H);
      
    % DM "nontradables" consumption share in total DM+CM consumption ... B
      
      qDM = ((SIGMA*(1-TAU_H))/A)...
                  *REL_PRICE*ZBAR*((K/H)^ALPHA)...
                    *(KAPPA*c_q(q,K,PHI)*q+(1-KAPPA)*c_q(q_czech,K,PHI)*q_czech);
              
      f6 = NTS - qDM/(X + qDM);
    
    % Collect equations
    lconval = [];                           % Linear binding constraints
    nonlconval = [f1,f2,f3,f6];      % Nonlinear binding constraints