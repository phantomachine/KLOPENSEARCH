function [mval, kval] = keulernb_bk(x,A,B,C,ETA,SIGMA,KAPPA,ALPHA,BETA,GAMA,...
                                 DELTA,ZBAR,EPSILON,THETA,THETA_B,b,...
                                    OMEGA_I,OMEGA_F,PHI,TAU_X,TAU_K,TAU_H)

    % KEULERNB.M
    % ---------------------------------------------------------------------
    % Function specifies the Euler equation in terms of steady state k:=K/H
    % for the Nash Bargaining case in DM model.
    %
    % Note: Compared with the PT and PROP cases, now we cannot write q as
    % an analytic function of K(k). So both (q, k) have to be determined
    % simultaneously using FSOLVE.
    %
    % Output : [qval, kval] = keulernb(...)
    %
    % ---------------------------------------------------------------------
    % (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
    %
    % See also NB_RUN, NB_EQ, NB_SS, SSMAPSTATIC_NASHBARG, FSOLVE,
    %          U_X, U_Q, G_NASH_Q, GAMMA_NASH, GAMMA_PROP
    % ---------------------------------------------------------------------
    
    k = x(1);
    q = x(2);
    
    % Composite parameters:
    
    REL_Pyh = (OMEGA_I^((EPSILON-1)/EPSILON))*THETA; % Ph/P
    
    % WRITE (q, q_czech, K, X) as functions of little k ...
    
    K = (((1+TAU_X)*(1-TAU_H)*B*REL_Pyh/A)*ZBAR*k^ALPHA)^(1/GAMA)...
          /(((1+TAU_X)^(-1))...
            *(OMEGA_I/OMEGA_F-((1-ALPHA)*TAU_H+ALPHA*TAU_K)*REL_Pyh)...
            *ZBAR*k^(ALPHA-1) - (1-TAU_K)*DELTA);
    
    q_czech = ((C*(ZBAR*K)^(PHI-1))/PHI)^(1/(ETA+PHI-1));
    
    X = ((1+TAU_X)*(1-TAU_H)*B*REL_Pyh*ZBAR*k^ALPHA /A)^(1/GAMA);
    
    % Now pin down q implicitly in steady-state Money Euler equation:
        uq = u_q(q,ETA,b,C);
        gq = g_nash_q(q,K,ZBAR,THETA_B,PHI,ETA,b,C);
        
    mval = - (1/SIGMA*KAPPA)*(1/BETA-1+SIGMA*KAPPA) * gq + uq;
    
    % Capital Euler equation at steady state to pin down k:
    
    gamma_NASH = gamma_nash(q,K,ZBAR,THETA_B,PHI,ETA,b,C);
    gamma_PROP = gamma_prop(q_czech,K,ZBAR,THETA_B,PHI,ETA,b,C);
    
    kval = - DELTA ...
            + (1-1/BETA)/(1-TAU_K) + REL_Pyh*ZBAR*ALPHA*k^(ALPHA-1) ...
                - ( SIGMA*(1+TAU_X)/((1-TAU_K)*U_X(X,B,GAMA))) ...
                  *( KAPPA*gamma_NASH + (1-KAPPA)*(1-THETA_B)*gamma_PROP );