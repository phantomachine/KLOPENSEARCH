function fun = ssk(kval,model)

% % SSK.m
% % =======================================================================
% % Price-Taking case. This function defines the Euler equation for optimal 
% % K at steady state, where it is re-written with a change of variables 
% % k := K/H. So ssk(k) is solely a univariate function in k.
% % Note that as:
% %         (i)  k --> 0, ssk(k) + 1 --> +infty
% %         (ii) k --> +infty, ssk(k) + 1 --> c, where 0 < c < 1.
% % So this implies that there exists a unique k* such that ssk(k*) = 0.
% % =======================================================================

    % % Extract variables 
        k = kval;
        
    % % Unpack parameters and steady state values needed
    	A = model.param.A;
        ETA1 = model.param.ETA1;
        ETA2 = model.param.ETA2;
        XI = model.param.XI;
        B = model.param.B;
        C = model.param.C;
        BETA = model.param.BETA;
        SIGMA = model.param.SIGMA;
        ZBAR = model.param.ZBAR;
        DELTA = model.param.DELTA;
        ALFA = model.param.ALFA;
        sx = model.ss.sx;
        sg = model.ss.sg;
  
    % % Define function:
    f1k = BETA*(1 + ALFA*ZBAR*k.^(ALFA-1) - DELTA);

        a0 = XI+ETA2-1;
        a1 = ALFA*a0;
        a2 = XI*ETA2*(ALFA-1);

    constant = BETA*(SIGMA*(XI-1)*(1-ALFA)*ZBAR/A)...
                          *(SIGMA*C/(XI*(1/BETA-1+SIGMA)))^(XI/a0);

    f2k = ( k.^((a1+a2)/a0) )...
                .*( (ZBAR - DELTA*k.^(1-ALFA))./((1+sg/sx)*((1/A)...
                    *(1-ALFA)*B*ZBAR*k.^ALFA).^(1/ETA1)) ).^(XI*ETA2/a0);

fun = (f1k + constant*f2k) - 1;  % ssk(k)