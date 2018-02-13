function model = pareto_eq(model)

% % PARETO_EQ.m
% % -----------------------------------------------------------------------
% % This function specifies equilibrium (functional) equations, transition
% % laws, and identities. In this case these are optimality conditions for
% % Pareto allocations. Note: 1 :=: Home and 2:=: Foreign.
% % 
% % INPUT: none
% % 
% % OUTPUT: model -- an OBJECT containing: 
% %             .x   : current states
% %             .xp  : next-period states
% %             .y   : current control
% %             .yp  : next-period control
% %             .f   : vector of equilibrium equations at zero
% %
% % REQUIRES: MATLAB Symbolic Math Toolbox
% %
% % -----------------------------------------------------------------------
% % Acknowledgement:
% % Based on algorithms by (c) Stephanie Schmitt-Grohe and Martin Uribe
% % -----------------------------------------------------------------------
% % (c), 2009 - , Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------

% % SYMBOLIC DEFINITIONS:

    % Define parameters
    syms   ALPHA BETA GAMA ETA DELTA EPSILON THETA PHI SIGMA ...
           A B b C ...
           ZBAR SG ...
           RHO_Z1   RHO_Z12...
           RHO_Z21  RHO_Z2 ...
         
    % Define variables, UPDATE: TOTX and TOTXp, 27/07/2010
    syms Z1 G1 X1 H1 q1 K1 yh1 yf1 Y1 I1        ... % Current Home
         Z2 G2 X2 H2 q2 K2 yf2 yh2 Y2 I2        ... % Current Foreign
         RELX TOTX RER GDP1 GDP2 CON1 CON2 INV1 INV2 ... % Current RER
         pq1 pq2 pX1 pX2 YPI1 YPI2 ...
         Z1p G1p X1p H1p q1p K1p yh1p yf1p Y1p I1p ... % Next-period Home
         Z2p G2p X2p H2p q2p K2p yf2p yh2p Y2p I2p ... % Next-period Foreign
         RELXp TOTXp RERp      ...                                    % Next-period RER
            GDP1p GDP2p CON1p CON2p INV1p INV2p ...
            pq1p pq2p pX1p pX2p YPI1p YPI2p
% % PRIMITIVE FUNCTIONS

    % Final goods aggregator
    GG1 = GG(yh1, yf1, EPSILON, THETA);
    GG2 = GG(yf2, yh2, EPSILON, THETA);
    
    % Intermediate goods production output
    FF1 = Z1*F(K1,H1,ALPHA);
    FF2 = Z2*F(K2,H2,ALPHA);
    
% % DERIVATIVE FUNCTIONS -  see directory: _func   

    % Marginal product of labor
    w1 = Z1*F_h(K1,H1,ALPHA);
    w2 = Z2*F_h(K2,H2,ALPHA);

    % Marginal utility of X
    MU_X1 = U_X(X1,B,GAMA);
    MU_X2 = U_X(X2,B,GAMA);
    
    % Marginal utility of q
    MU_q1 = u_q(q1,ETA,b,C);
    MU_q2 = u_q(q2,ETA,b,C);
    
    % Marginal cost of q
    MC_q1 = c_q(q1,K1,PHI) / (Z1^PHI); % Added Z1: 17/8/2010
    MC_q2 = c_q(q2,K2,PHI) / (Z2^PHI); % Added Z2: 17/8/2010
    
    
    % Marginal products wrt yh1 and yf1
    GG_yh1 = GG_1(yh1,yf1,EPSILON,THETA);
    GG_yf1 = GG_2(yh1,yf1,EPSILON,THETA);
    
    % Marginal products wrt yf2 and yh2
    GG_yf2 = GG_1(yf2,yh2,EPSILON,THETA);
    GG_yh2 = GG_2(yf2,yh2,EPSILON,THETA);
    
% % ONE-PERIOD-AHEAD R.V.'s ("p" suffix for "+")
    
    % Marginal product of capital "+"
    r1p = Z1p*F_k(K1p,H1p,ALPHA);
    r2p = Z2p*F_k(K2p,H2p,ALPHA);
    
    % Marginal utility of X "+"
    MU_X1p = U_X(X1p,B,GAMA);
    MU_X2p = U_X(X2p,B,GAMA);
    
    % Marginal cost wrt K "+"
    MC_K1p = c_k(q1p,K1p,PHI) / (Z1p^PHI); % Added Z1p: 17/8/2010;
    MC_K2p = c_k(q2p,K2p,PHI) / (Z2p^PHI); % Added Z2p: 17/8/2010;
    
    % Marginal products wrt yh1 "+"
    GG_yh1p = GG_1(yh1p,yf1p,EPSILON,THETA);
    
    % Marginal products wrt yf2 "+"
    GG_yf2p = GG_1(yf2p,yh2p,EPSILON,THETA);
    
% % DEFINE EQUILIBRIUM/OPTIMALITY FUNCTIONAL EQUATIONS HERE:

    % Country 1 (Home) Euler Equation
    f1 = MU_X1 - BETA*( MU_X1p*(GG_yh1p*r1p + 1 - DELTA) ...
                                                - SIGMA*MC_K1p);
    f2 = MU_X2 - BETA*( MU_X2p*(GG_yf2p*r2p + 1 - DELTA) ...
                                                - SIGMA*MC_K2p);
    % Efficient q allocation
    f3 = MU_q1 - MC_q1;
    f4 = MU_q2 - MC_q2;
    
    % Real exchange rate or International CM Terms of Trade
    f5 = TOTX - GG_yf1/GG_yh1;                  % UPDATE: TOTX, 27/07/2010
    
    % Efficient H allocation
    f6 = A - w1*MU_X1*GG_yh1;
    f7 = A - w2*MU_X2*GG_yf2;
    
    % Shadow price of yh and yf
    f8 = GG_yh1*MU_X1 - GG_yh2*MU_X2;
    f9 = GG_yf1*MU_X1 - GG_yf2*MU_X2;
    
    % Feasibility constraint for yh and yf
    f10 = FF1 - yh1 - yh2;
    f11 = FF2 - yf2 - yf1;
    
    % Resource constraint
    f12 = GG1 - X1 - K1p + (1-DELTA)*K1 - G1;
    f13 = GG2 - X2 - K2p + (1-DELTA)*K2 - G2;
    
    % Identitities
    
    pq1 = MU_X2 / MU_q1; % New: 15/08/2010
    pq2 = MU_X2 / MU_q2;
    
    pX1 = TOTX; % TOTX = MU_X2 / MU_X1;
    pX2 = 1;    % numeraire is X2
        
    f14 = Y1 - ((yh1 + yh2)*GG_yh1 - SIGMA*pq1*q1)/pX1; %Y1 - GG1 - SIGMA*q1; % ERROR:22-02-10
    f15 = Y2 - ((yf2 + yf1)*GG_yf2 - SIGMA*pq2*q2)/pX2; %Y2 - GG2 - SIGMA*q2; % ERROR:22-02-10
    
    % Non-traded (q) good share in total consumption
    NTS1 = SIGMA*pq1*q1 / (SIGMA*pq1*q1 + pX1*Y1);
    NTS2 = SIGMA*pq2*q2 / (SIGMA*pq2*q2 + pX2*Y2);
    
    % Investment
    f16 = I1 - K1p + (1-DELTA)*K1;
    f17 = I2 - K2p + (1-DELTA)*K2;
    
    % Relative CM consumption
    f18 = RELX - CON1/CON2;
    
    
    % RER definition in non-monetary economy/Pareto allocation % UPDATE:
    % TOTX, 27/07/2010
            
    f28 = YPI1 - NTS1*pq1 - (1-NTS1)*pX1;   % GDP deflator : 15/08/2010
    f29 = YPI2 - NTS2*pq2 - (1-NTS2)*pX2;
    
    f19 = RER - YPI2/YPI1;
    
    f20 = GDP1 - Y1*pX1/YPI1;
   
    f21 = GDP2 - Y2*pX2/YPI2;
    
    f22 = CON1 - (pX1*X1 + SIGMA*pq1*q1)/YPI1;
    
    f23 = CON2 - (pX2*X2 + SIGMA*pq2*q2)/YPI2; 
    
    f24 = INV1 - pX1*I1/YPI1;
    
    f25 = INV2 - pX2*I2/YPI2;
    
    
    % Exogenous shock processes
    f26 = log(Z1p) - RHO_Z1 * log(Z1) - RHO_Z12 * log(Z2); 
    f27 = log(Z2p) - RHO_Z2 * log(Z2) - RHO_Z21 * log(Z1);
    
  
   
    
% % STACK SYSTEM AS f
f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;...
                                                    f17;f18;f19;f20;...
                                                    f21;f22;f23;f24;f25;...
                                                    f28;f29;f26;f27];

    % Define the vector of states, x, and controls, y: 
    x = [ K1 K2 Z1 Z2 ];
    y = [ X1 X2 H1 H2 q1 q2 yh1 yf1 yf2 yh2 Y1 Y2 I1 I2 ...
                                TOTX RER RELX ... % UPDATE: TOTX, 27/07/2010
                                GDP1 GDP2 CON1 CON2 INV1 INV2 YPI1 YPI2 ];

    xp = [ K1p K2p Z1p Z2p ];
    yp = [ X1p X2p H1p H2p q1p q2p yh1p yf1p yf2p yh2p Y1p Y2p I1p I2p ...
                            TOTXp RERp RELXp  ... % UPDATE: TOTX, 27/07/2010
                                GDP1p GDP2p CON1p CON2p INV1p INV2p YPI1p YPI2p ];
% % Pack symbolic vectors into MODEL
% 
% model.parasym = [A B b C BETA DELTA ALFA THETA ETA1 ETA2 XI SIGMA ...
%                                       ZBAR GBAR ...
%                                       RHO_Z RHO_G ...
%                                       RHO_Zstar RHO_Gstar ...
%                                       RHO_ZZstar  RHO_ZstarZ ...
%                                       RHO_PSI RHO_PSIstar ];
model.x = x;
model.y = y;
model.xp = xp;
model.yp = yp;
model.f = f; 