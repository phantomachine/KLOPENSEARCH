function model = pt_eq

% % PT_EQ.m
% % -----------------------------------------------------------------------
% % This function specifies equilibrium (functional) equations, transition
% % laws, and identities. In this case these are optimality conditions for
% % RCE allocations. Note: 1 :=: Home and 2:=: Foreign.
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
syms   ALPHA BETA GAMA ETA DELTA EPSILON THETA PHI SIGMA KAPPA ...
    TAU_X TAU_H TAU_K ...
    A B b C ...
    ZBAR SG ...
    RHO_Z1   RHO_Z12...
    RHO_Z21  RHO_Z2 ...
    MU_M1 MU_M12 ...                % Extra syms for monetary shocks
    MU_M21 MU_M2

% Define variables
syms Z1 psi1 G1 X1 H1 q1 qz1 K1 yh1 yf1 Y1 I1 NX1 ... % Current Home
    Z2 psi2 G2 X2 H2 q2 qz2 K2 yf2 yh2 Y2 I2 NX2 ... % Current Foreign
    Ptilde1 Ptilde2...                               % Prices of q's
    P1 Pyh1 Pyh2 ...                                 % Price levels
    P2 Pyf2 Pyf1 ...                                 % Price levels
    Ptilde1 Ptilde2...                               % DM prices
    YPI1 YPI2...                                     % ad-hock Price levels
    RER E RELX TOT RER_YPI NTS_s1 NTS_s2...                        % Current RER
    Agg_C1 Agg_C2...
    Agg_Y1 Agg_Y2...
    Agg_I1 Agg_I2 P2P1 ...
    Z1p psi1p G1p X1p H1p q1p qz1p K1p ...
    yh1p yf1p Y1p I1p NX1p...
    Z2p psi2p G2p X2p H2p q2p qz2p K2p ...
    yf2p yh2p Y2p I2p NX2p ...
    P1p Pyh1p Pyh2p ...                              % Price levels
    P2p Pyf2p Pyf1p ...                              % Price levels
    Ptilde1p Ptilde2p...                             % DM prices
    YPI1p YPI2p...                                     % ad-hock Price levels
    RERp Ep RELXp TOTp RER_YPIp...                   % Next-period RER
    Agg_C1p Agg_C2p...
    Agg_Y1p Agg_Y2p...
    Agg_I1p Agg_I2p P2P1p

% % PRIMITIVE FUNCTIONS

% Final goods aggregator
GG1 = GG(yh1, yf1, EPSILON, THETA);
GG2 = GG(yf2, yh2, EPSILON, THETA);

% Intermediate goods production output
FF1 = Z1*F(K1,H1,ALPHA);
FF2 = Z2*F(K2,H2,ALPHA);

% % DERIVATIVE FUNCTIONS -  see directory: _func

% Marginal products wrt yh1 and yf1
GG_yh1 = GG_1(yh1,yf1,EPSILON,THETA);
%GG_yf1 = GG_2(yh1,yf1,EPSILON,THETA);

% Marginal products wrt yf2 and yh2
GG_yf2 = GG_1(yf2,yh2,EPSILON,THETA);
%GG_yh2 = GG_2(yf2,yh2,EPSILON,THETA);

% Rental:  labor
w1 = GG_yh1*Z1*F_h(K1,H1,ALPHA);
w2 = GG_yf2*Z2*F_h(K2,H2,ALPHA);

% Rental:  capital "+"
r1 = GG_yh1*Z1*F_k(K1,H1,ALPHA);
r2 = GG_yf2*Z2*F_k(K2,H2,ALPHA);

% Marginal utility of X
MU_X1 = U_X(X1,B,GAMA);
MU_X2 = U_X(X2,B,GAMA);

% % ONE-PERIOD-AHEAD R.V.'s ("p" suffix for "+")

% Marginal products wrt yh1 "+"
GG_yh1p = GG_1(yh1p,yf1p,EPSILON,THETA);

% Marginal products wrt yf2 "+"
GG_yf2p = GG_1(yf2p,yh2p,EPSILON,THETA);

% Marginal product of capital "+"
r1p = GG_yh1p*Z1p*F_k(K1p,H1p,ALPHA);
r2p = GG_yf2p*Z2p*F_k(K2p,H2p,ALPHA);

% Marginal utility of X "+"
MU_X1p = U_X(X1p,B,GAMA);
MU_X2p = U_X(X2p,B,GAMA);

% Marginal cost wrt K "+"
MC_K1p = c_k(q1p,K1p,PHI) / (Z1p^PHI); % Added Z1p .. 17/08/2010
MC_K2p = c_k(q2p,K2p,PHI) / (Z2p^PHI); % Added Z2p .. 17/08/2010

MC_K1zp = c_k(qz1p,K1p,PHI) / (Z1p^PHI); % Added Z1p .. 17/08/2010
MC_K2zp = c_k(qz2p,K2p,PHI) / (Z1p^PHI); % Added Z2p .. 17/08/2010

% Marginal utility of q "+"
MU_q1p = u_q(q1p,ETA,b,C);
MU_q2p = u_q(q2p,ETA,b,C);

% Marginal cost of q "+"
MC_q1p = c_q(q1p,K1p,PHI) / (Z1p^PHI); % Added Z1 .. 17/08/2010;
MC_q2p = c_q(q2p,K2p,PHI) / (Z2p^PHI); % Added Z2 .. 17/08/2010;

% Marginal cost of q "+"
MC_q1 = c_q(q1,K1,PHI) / (Z1^PHI); % Added Z1 .. 17/08/2010;
MC_q2 = c_q(q2,K2,PHI) / (Z2^PHI); % Added Z2 .. 17/08/2010;

% Marginal utility of qz (q_czech)     % Added 04/05/2010
MU_qz1 = u_q(qz1,ETA,b,C);
MU_qz2 = u_q(qz2,ETA,b,C);

% Marginal cost of qz (q_czech)        % Added 04/05/2010
MC_qz1 = c_q(qz1,K1,PHI) / (Z1^PHI); % Added Z1p .. 17/08/2010
MC_qz2 = c_q(qz2,K2,PHI) / (Z2^PHI); % Added Z1p .. 17/08/2010


% % DEFINE EQUILIBRIUM/OPTIMALITY FUNCTIONAL EQUATIONS HERE:

% Note: Functions are stored in directory /_func

% Country 1-2 (Home-Foreign) Euler Equation, FOC K1+, K2+

f1 = MU_X1 - BETA*( MU_X1p*( 1 + (1-TAU_K)*(r1p - DELTA) ) ...
    - SIGMA*(1+TAU_X)*(KAPPA*MC_K1p + (1-KAPPA)*MC_K1zp) );
f2 = MU_X2 - BETA*( MU_X2p*( 1 + (1-TAU_K)*(r2p - DELTA) ) ...
    - SIGMA*(1+TAU_X)*(KAPPA*MC_K2p + (1-KAPPA)*MC_K2zp) );
% Pins down (K1+, K2+)
% KAPPA: Added 04/05/2010

% Country 1-2 (Home-Foreign), Euler equation, FOC, M1+, M2+

%     f3 = g(q1,K1,PHI) ...
%         - BETA*g(q1p,K1p,PHI)*(1/psi1)*( 1 - SIGMA + SIGMA*MU_q1p/MC_q1p );
%     f4 = g(q2,K2,PHI) ...
%         - BETA*g(q1p,K1p,PHI)*(1/psi2)*( 1 - SIGMA + SIGMA*MU_q2p/MC_q2p );

f3 = MU_X1 ...
    - BETA*MU_X1p*(P1/P1p)*(1/psi1)...
    *( 1 - SIGMA*KAPPA + SIGMA*KAPPA*MU_q1p/MC_q1p );
f4 = MU_X2 ...
    - BETA*MU_X2p*(P2/P2p)*(1/psi2)...
    *( 1 - SIGMA*KAPPA + SIGMA*KAPPA*MU_q2p/MC_q2p );
% Pins down (q1, q2)
% KAPPA: Added 04/05/2010

% Country 1-2 (Home-Foreign) Labor market clearing

f5 = A*(1+TAU_X) - (1-TAU_H)*w1*MU_X1;
f6 = A*(1+TAU_X) - (1-TAU_H)*w2*MU_X2;       % Pins down (X1, X2)

% Demand for intermediate goods (yh1, yh2)

f7 = yh1 - demand(P1,Pyh1,GG1,EPSILON,THETA);
f8 = yh2 - demand(P2,Pyh2,GG2,EPSILON,1-THETA); % Pins down (yh1, yh2)

% Demand for intermediate goods (yh1, yh2)
f9 = yf1 - demand(P1,Pyf1,GG1,EPSILON,1-THETA);
f10 = yf2 - demand(P2,Pyf2,GG2,EPSILON,THETA);  % Pins down (yf1, yf2)

% Market clearing for yh and yf
f11 = FF1 - yh1 - yh2;
f12 = FF2 - yf2 - yf1;                          % Pins down (Ph1, Pf1)

% Resource constraint
G1 = (TAU_X*X1 + TAU_H*w1*H1 + TAU_K*(r1-DELTA)*K1);
G2 = (TAU_X*X2 + TAU_H*w2*H2 + TAU_K*(r2-DELTA)*K2);

f13 = GG1 - X1 - K1p + (1-DELTA)*K1 - G1;
f14 = GG2 - X2 - K2p + (1-DELTA)*K2 - G2;       % Pins down (X1, X2)

% Identitities - output GDP
%     f15 = Y1 - GG1 - Pyh2*E*yh2/P1 + Pyf1*yf1/P1 - SIGMA*KAPPA*q1...
%                                                 - SIGMA*(1-KAPPA)*qz1;
%     f16 = Y2 - GG2 - Pyf1*yf1/P2 + Pyh2*E*yh2/P2 - SIGMA*KAPPA*q2...
%                                                 - SIGMA*(1-KAPPA)*qz2;
Ptilde_q1 = ((A*psi1/(P1*(1-TAU_H)*w1))^(-1))*MC_q1;
% Ptilde_q1 = c_q(q1,K1,PHI);
% Home money trades: price
% NOTE: c_q, Marginal Cost in utils
Ptilde_qz1 = ((A*psi1/(P1*(1-TAU_H)*w1))^(-1))*MC_qz1;
% Ptilde_qz1 = c_q(qz1,K1,PHI);
% Home credit trades: price
Ptilde_q2 = ((A*psi2/(P2*(1-TAU_H)*w2))^(-1))*MC_q2;
% Ptilde_q2 = c_q(q2,K2,PHI);
% Foreign money trades: price
Ptilde_qz2 = ((A*psi2/(P2*(1-TAU_H)*w2))^(-1))*MC_qz2;
% Ptilde_qz2 = c_q(qz2,K2,PHI);
% Foreign credit trades: price

% Define the share of DM output
NTS_s1 = 1-Pyh1*FF1/(Pyh1*FF1+SIGMA*(KAPPA*Ptilde_q1*q1+(1-KAPPA)*Ptilde_qz1*qz1));
NTS_s2 = 1-Pyf2*FF2/(Pyh2*FF2+SIGMA*(KAPPA*Ptilde_q2*q2+(1-KAPPA)*Ptilde_qz2*qz2));

% f15 = Y1 - (Pyh1/P1)*FF1 - SIGMA*KAPPA*q1...
%     - SIGMA*(1-KAPPA)*qz1;
% f16 = Y2 - (Pyf2/P2)*FF2 - SIGMA*KAPPA*q2...
%     - SIGMA*(1-KAPPA)*qz2;
f15 = P1*Y1 - Pyh1*FF1 - SIGMA*KAPPA*Ptilde_q1*q1...
    - SIGMA*(1-KAPPA)*Ptilde_qz1*qz1;
f16 = P2*Y2 - Pyf2*FF2 - SIGMA*KAPPA*Ptilde_q2*q2...
    - SIGMA*(1-KAPPA)*Ptilde_qz2*qz2;
% Pins down (Y1, Y2)

%     % Identitities - CPI level
%     f17 = P1 - cpi(Ph1,Pf1,EPSILON,THETA);
%     f18 = P2 - cpi(Pf2,Ph2,EPSILON,THETA);

% DM FOC pricing outcome
f17 = A*psi1/(P1*(1-TAU_H)*w1) - g(q1, K1, Z1, PHI); % Added Z1.. 17/08/2010
f18 = A*psi2/(P2*(1-TAU_H)*w2) - g(q2, K2, Z2, PHI); % Added Z2.. 17/08/2010
% Pins down (P1, P2)
% DM credit exhcange - First best allocation
f19 = MU_qz1 - MC_qz1;
f20 = MU_qz2 - MC_qz2;                           % Pins down (qz1, qz2)
% Added 04/05/2010

% Identitities - Law of one price
f21 = Pyh1 - E*Pyh2;
f22 = Pyf1 - E*Pyf2;                             % Pins down (Ph2, Pf2)

% Identitities - Investment
f23 = I1 - K1p + (1-DELTA)*K1;
f24 = I2 - K2p + (1-DELTA)*K2;                   % Pins down (I1, I2)


% Real exchange rate or International CM Terms of Trade
f25 = RER - E * YPI2 / YPI1;                        % Pins down RER

% Nominal ER -  risk sharing complete markets
% C.f.: f24 = (MU_X1p/MU_X1)*(P1/P1p) ...
%                - (MU_X2p/MU_X2)*(P2/P2p)*(E/Ep); % incomplete markets

f26 = E*P2/P1 - MU_X2/MU_X1;              % complete markets (CM bonds)
% Pins down E

% Relative CM consumption
f27 = RELX - Agg_C1/Agg_C2;

% Identitities - Net export
f28 = NX1 - Pyh2*E*yh2/P1 + Pyf1*yf1/P1;
f29 = NX2 - Pyf1*yf1/P2 + Pyh2*E*yh2/P2;        % Pins down (NX1, NX2)

f34 = TOT - Pyf1/Pyh1;

%     % Constructing CPI  added 27/05/2010
%     f35 = Ptilde1 - c_q(qz1,K1,PHI);
%     f36 = Ptilde2 - c_q(qz2,K1,PHI);
%     f37 = CPI1 - ((X1/(X1+SIGMA*qz1))*P1+(SIGMA*qz1/(X1+SIGMA*qz1))*Ptilde1);
%     f38 = CPI2 - ((X2/(X2+SIGMA*qz2))*P2+(SIGMA*qz2/(X2+SIGMA*qz2))*Ptilde2);

% Constructing CPI  added 27/05/2010

f35 = Ptilde1 - ( KAPPA*Ptilde_q1 + (1-KAPPA)*Ptilde_qz1 );
% Average DM price level: Home

f36 = Ptilde2 - ( KAPPA*Ptilde_q2 + (1-KAPPA)*Ptilde_qz2 );
% Average DM price level: Foreign

% qavg1 = KAPPA*q1 + (1-KAPPA)*qz1;
f37 = YPI1 - (1-NTS_s1)*P1 - NTS_s1*Ptilde1;

% qavg2 = KAPPA*q2 + (1-KAPPA)*qz2;
f38 = YPI2 - (1-NTS_s2)*P2 - NTS_s2*Ptilde2;


f39 = RER_YPI - E* YPI2 / YPI1;

% define Aggregate consumption
f40 = Agg_C1 - ((P1*X1+SIGMA*(1-KAPPA)*Ptilde_qz1*qz1+...
                          SIGMA*KAPPA*Ptilde_q1*q1)/YPI1);
f41 = Agg_C2 - ((P2*X2+SIGMA*(1-KAPPA)*Ptilde_qz2*qz2+...
                          SIGMA*KAPPA*Ptilde_q2*q2)/YPI2);
f42 = Agg_Y1 - (Pyh1*FF1+SIGMA*(KAPPA*Ptilde_q1*q1+(1-KAPPA)*Ptilde_qz1*qz1))/YPI1;
f43 = Agg_Y2 - (Pyf2*FF2+SIGMA*(KAPPA*Ptilde_q2*q2+(1-KAPPA)*Ptilde_qz2*qz2))/YPI2;

f44 = Agg_I1 - P1*I1/YPI1;
f45 = Agg_I2 - P2*I2/YPI2;

f46 = P2P1 - YPI2/YPI1;

% Exogenous shock processes
f30 = log(Z1p) - RHO_Z1 * log(Z1) - RHO_Z12 * log(Z2);
f31 = log(Z2p) - RHO_Z2 * log(Z2) - RHO_Z21 * log(Z1);
f32 = log(psi1p) - MU_M1 * log(psi1) - MU_M12 * log(psi2);
f33 = log(psi2p) - MU_M2 * log(psi2) - MU_M21 * log(psi1);


% % STACK SYSTEM AS f
f = [   f1;f2;f3;f4;...
    f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;f17;f18;...
    f19;f20;f21;f22;f23; f24;f25;f26; ...
    f27;f28; f29; f30; f31; f32; f33; f34 ; ...
    f35;f36;f37;f38;f39;f40;f41;f42;f43;f44;f45; f46  ];

% Define the vector of states, x, and controls, y:
x = [ K1 K2 Z1 Z2 psi1 psi2];
y = [ X1 X2 H1 H2 q1 q2 qz1 qz2 yh1 yf1 yf2 yh2 Y1 Y2 ...
    E I1 I2 P1 Pyh1 Pyh2 P2 Pyf2 Pyf1 RER RELX NX1 NX2 TOT...
    Ptilde1 Ptilde2 YPI1 YPI2 RER_YPI Agg_C1 Agg_C2 Agg_Y1 Agg_Y2...
     Agg_I1 Agg_I2 P2P1];

xp = [ K1p K2p Z1p Z2p psi1p psi2p];
yp = [ X1p X2p H1p H2p q1p q2p qz1p qz2p yh1p yf1p yf2p yh2p Y1p Y2p ...
    Ep I1p I2p P1p Pyh1p Pyh2p P2p Pyf2p Pyf1p RERp RELXp NX1p NX2p TOTp...
    Ptilde1p Ptilde2p YPI1p YPI2p RER_YPIp Agg_C1p Agg_C2p Agg_Y1p Agg_Y2p...
    Agg_I1p Agg_I2p P2P1p];

% % Pack symbolic vectors into MODEL

model.x = x;
model.y = y;
model.xp = xp;
model.yp = yp;
model.f = f;