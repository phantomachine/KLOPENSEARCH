function fval = keulerpt(k,A,B,C,ETA,SIGMA,KAPPA,ALPHA,BETA,GAMA,DELTA,ZBAR,EPSILON,THETA,...
                            OMEGA_I,OMEGA_F,PHI,TAU_X,TAU_K,TAU_H)

    % KEULERPT.M
    % Function specifies the Euler equation in terms of steady state k:=K/H

    % Composite parameters:
    % ---------------------------------------------------------------------
    
    REL_Pyh = (OMEGA_I^((EPSILON-1)/EPSILON))*THETA; % Ph/P
    
    % WRITE (q, q_czech, K, X) as functions of little k ...
    
    K = (((1+TAU_X)*(1-TAU_H)*B*REL_Pyh/A)*ZBAR*k^ALPHA)^(1/GAMA)...
          /(((1+TAU_X)^(-1))...
            *(OMEGA_I/OMEGA_F-((1-ALPHA)*TAU_H+ALPHA*TAU_K)*REL_Pyh)...
            *ZBAR*k^(ALPHA-1) - (1-TAU_K)*DELTA);
    
    q_czech = ((C*K^(PHI-1))/PHI)^(1/(ETA+PHI-1));
    
    q = ((SIGMA*KAPPA*C*K^(PHI-1))...
                    /(PHI*(1/BETA-1+SIGMA*KAPPA)) )^(1/(ETA+PHI-1));
                
    X = ((1+TAU_X)*(1-TAU_H)*B*REL_Pyh*ZBAR*k^ALPHA /A)^(1/GAMA);
    
    % Capital Euler equation at steady state
    
    fval = - DELTA ...
            + (1-1/BETA)/(1-TAU_K) + REL_Pyh*ZBAR*ALPHA*k^(ALPHA-1) ...
                - ( SIGMA*(1+TAU_X)/((1-TAU_K)*U_X(X,B,GAMA))) ...
                    *( KAPPA*c_k(q,K,PHI) + (1-KAPPA)*c_k(q_czech,K,PHI) );