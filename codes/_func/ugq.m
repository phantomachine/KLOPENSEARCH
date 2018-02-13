function [f1,f2] = ugq(qk,model)
% function [f1,f2] = ugq(q, KDM, model)
% This function specifies the Euler equation for optimal money holdings and
% capital which is a pair of functional equations in (q, KDM), evaluated at
% steady state.
% We need to solve for (q,KDM) simultaneously.
% =========================================================================
% (c) 2009, T.Kam and J.Lee

    % % Extract variables
        q = qk(1); 
        K = qk(2);
        
    % % Unpack parameters and steady state values needed
        %THETA = model.param.THETA;
        ETA1 = model.param.ETA1;
        ETA2 = model.param.ETA2;
        XI = model.param.XI;
        B = model.param.B;
        b = model.param.b;
        C = model.param.C;
        BETA = model.param.BETA;
        SIGMA = model.param.SIGMA;
        ZETABAR = model.param.ZETABAR;
        DELTA = model.param.DELTA;
        ZETA = ZETABAR;
        ALFA = model.param.ALFA;
        
        KYss = model.ss.KYss;
     
        X = model.ss.X;
 
    % % USER: Hard-code the primtive functions to save time
    f1 = (1/BETA - (1-SIGMA))*ZETA*c_q(q,K,XI) - SIGMA*u_q(q, ETA2, b, C);
    
    f2 = (1/BETA - (1 + ALFA/KYss - DELTA))*U_X(X,B,ETA1) + SIGMA*ZETA*c_k(q,K,XI);