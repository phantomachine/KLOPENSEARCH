function fval = keulercomp(k,ALPHA,BETA,GAMMA,DELTA,ZBAR,...
                            OMEGA_I,OMEGA_F,LAM_1,LAM_2,...
                                        LAM_3,LAM_4,LAM_5,LAM_6,LAM_7,CHI,SG)

    % KMAP.M
    % Function specifies the Euler equation in terms of steady state k:=K/H
   
    
    GHratio = SG*(OMEGA_I/OMEGA_F)*ZBAR*k.^ALPHA;
    
    RHS = LAM_1*k.^(ALPHA-1) ...
            + (LAM_2*k.^ALPHA).*( (LAM_3*k.^(ALPHA/GAMMA) + LAM_6*k.^ALPHA + LAM_7*k.^(ALPHA-1) )...
             ./ (LAM_4*k.^(ALPHA-1) - LAM_5*(1 + GHratio/DELTA)) ).^(-CHI);

    LHS = 1/BETA - 1 + DELTA;

    fval = RHS - LHS;