%function model = pareto_ss(model,SWEEP1,varargin)

% % PARETO_SS.m
% % -----------------------------------------------------------------------
% % This function compute the planner's steady state allocations.
% % 
% % INPUT: optional
% % 
% % OUTPUT: model.ss (containing parameters and steady-state values)
% %
% % REMARK: This is also where the USER defines the deep parameter values.
% %
% % -----------------------------------------------------------------------
% % Acknowledgement:
% % Based of algorithms by (c) Stephanie Schmitt-Grohe and Martin Uribe
% % -----------------------------------------------------------------------
% % (c) 2009 - , Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------

% % =======================================================================
% % STEP 1: Set numeric values for deep parameters
% %         REMARK: Taken from "Model 2" calibration in AWW Table 2. Some
% %                 parameters are normalized to a quarterly frequency.
% % =======================================================================

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('PARETO_SS.M: STEP 1. SET PARAMETERS ...')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')

% INVOKE PARAMETER SETTINGS SCRIPT:
parameterscommon;

% Calibration from pt_ss ...
% (Pricetaking baseline with credit, money and DM capital):
        ALPHA = 0.33;
        A = 0.3326;
        B = 0.0967;
        SIGMA = 0.078;
        PHI = 1.289;   
        %PHI = 1.000000000001;


LAM_1 = (OMEGA_I^((EPSILON-1)/EPSILON)) * THETA * ALPHA * ZBAR;

LAM_2 = (1/A)*( SIGMA*(PHI-1)*(C/PHI)^(PHI/(ETA+PHI-1)) ) ...
                    *(1-ALPHA)*ZBAR*THETA*OMEGA_I^((EPSILON-1)/EPSILON);

LAM_3 = (ZBAR*(1-ALPHA)*B*THETA*OMEGA_I^((EPSILON-1)/EPSILON) )^(1/GAMA);

LAM_4 = (OMEGA_I/OMEGA_F)*ZBAR*A^(1/GAMA);

LAM_5 = DELTA*A^(1/GAMA);

CHI = PHI*ETA/(ETA + PHI - 1);

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('PARETO_SS.M: STEP 2. SOLVING FOR PLANNER STEADY-STATE ALLOCATION')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')

% Solve Euler Equation for k := K/H
options0 = optimset( 'Display','off',...
                     'TolFun',10e-12,...
                     'MaxFunEvals',10e+5,...
                     'MaxIter',10e+5         );

k_init = 0.1;
                
[kss, fss, flag, model.fsout] = fsolve('keuler',k_init,options0,...
                                    ALPHA,BETA,GAMA,DELTA,ZBAR,...
                                    OMEGA_I,OMEGA_F,LAM_1,LAM_2,...
                                                LAM_3,LAM_4,LAM_5,CHI,SG);

if flag && strcmp(options0.Display,'off')   
    fprintf('FSOLVE: Equation Solved.\n\t Type model.fsout for more...\n');  
else 
    warning('Error:FSOLVE','You may not have true zero')
end

if strcmp(options0.Display,'on')
    % Check visually
    k = linspace(0,100,100);
    fval = keuler(k,ALPHA,BETA,GAMA,DELTA,ZBAR,...
                            OMEGA_I,OMEGA_F,LAM_1,LAM_2,...
                                        LAM_3,LAM_4,LAM_5,CHI,SG);
    figure
        hold on
        plot(k,fval,k,zeros(size(k)));
        plot(kss, fss, 'o','MarkerFaceColor','r');
        text(kss,fss, sprintf('(%0.2g,%0.2g)',kss,fss) );
        hold off
        xlabel('k')
end

% Set 1 :=: Home, 2 :=: Foreign
kss1 = kss;

% Solve for Xss as function of kss
Xss1 = (ZBAR*(1-ALPHA)*B*THETA*(OMEGA_I^((EPSILON-1)/EPSILON))/A)...
                                                               *kss1^ALPHA;

% Solve Kss as function of kss
GHratio = SG*(OMEGA_I/OMEGA_F)*ZBAR*kss1^ALPHA;
Kss1 = Xss1 /((OMEGA_I/OMEGA_F)*ZBAR*kss1^(ALPHA-1) -DELTA - GHratio);

% Solve for Hss
Hss1 = Kss1/kss1;

% Solve for Gss
Gss1 = Hss1*GHratio;

% Solve for qss as function of Kss
qss1 = ((C/PHI)^(1/(ETA+PHI-1)))*(Kss1^((PHI-1)/(ETA+PHI-1)));

% Define output of Home intermediate good sector
Yss1 = ZBAR*(Kss1^(ALPHA))*(Hss1^(1-ALPHA));

% Solve for yhss1
yhss1 = Yss1/OMEGA_F;

% Solve for yfss1
yfss1 = yhss1*(THETA/(1-THETA))^(EPSILON/(1-EPSILON));

% Iss1
Iss1 = DELTA*Kss1;

% By symmetry, for Foreign's Pareto allocation:
kss2 = kss1;
Xss2 = Xss1;
Kss2 = Kss1;
Hss2 = Hss1;
qss2 = qss1;
Yss2 = Yss1;
yhss2 = yfss1;
yfss2 = yhss1;
Gss2 = Gss1;
Zss1 = ZBAR;
Zss2 = Zss1;
Iss2 = Iss1;

% Marginal products wrt yh1 and yf1
    GG_yh1 = GG_1(yhss1,yfss1,EPSILON,THETA);
    GG_yf1 = GG_2(yhss1,yfss1,EPSILON,THETA);
    
    % Marginal products wrt yf2 and yh2
    GG_yf2 = GG_1(yfss2,yhss2,EPSILON,THETA);
    GG_yh2 = GG_2(yfss2,yhss2,EPSILON,THETA);


% Steady state real exchange rate
TOTX = 1;
%RER = 1;
RELX = 1;

 % Marginal utility of X
    MU_X1 = U_X(Xss1,B,GAMA);
    MU_X2 = U_X(Xss2,B,GAMA);
    
    % Marginal utility of q
    MU_q1 = u_q(qss1,ETA,b,C);
    MU_q2 = u_q(qss2,ETA,b,C);

    pq1 = MU_X2 / MU_q1; % New: 15/08/2010
    pq2 = MU_X2 / MU_q2;
    
    pX1 = TOTX; % TOTX = MU_X2 / MU_X1;
    pX2 = 1;    % numeraire is X2
        
    Y1 = (yhss1 + yhss2)*GG_yh1 - SIGMA*pq1*qss1/pX1; %Y1 - GG1 - SIGMA*q1; % ERROR:22-02-10
    Y2 = (yfss2 + yfss1)*GG_yf2 - SIGMA*pq2*qss2/pX2; %Y2 - GG2 - SIGMA*q2; % ERROR:22-02-10

    % Non-traded (q) good share in total consumption
    NTS1 = SIGMA*pq1*qss1 / (SIGMA*pq1*qss1 + pX1*Y1);
    NTS2 = SIGMA*pq2*qss2 / (SIGMA*pq2*qss2 + pX2*Y2);
    
             
        YPI1 = NTS1*pq1 + (1-NTS1)*pX1;   % GDP deflator : 15/08/2010
        YPI2 = NTS2*pq2 + (1-NTS2)*pX2;
    
    RER = YPI2/YPI1;
    
    GDP1 = Y1*pX1/YPI1;
   
    GDP2 = Y2*pX2/YPI2;
    
    CON1 = (pX1*Xss1 + SIGMA*pq1*qss1)/YPI1;
    
    CON2 = (pX2*Xss2 + SIGMA*pq2*qss2)/YPI2; 
    
    INV1 = pX1*Iss1/YPI1;
    
    INV2 = pX2*Iss2/YPI2;

% STORE PARAMETERS IN OBJECT
model.params = {     ALPHA      ... % 1
                     BETA       ... % 2
                     GAMA      ... % 3
                     DELTA      ... % 4
                     PHI        ... % 5
                     ETA        ... % 6
                     SIGMA      ... % 7
                     THETA      ... % 8
                     EPSILON    ... % 9
                     ZBAR       ... % 10
                     A          ... % 11
                     B          ... % 12
                     C          ... % 13
                     b          ... % 14
                     RHO        ... % 15
                     ECOV       ... % 16
                     OMEGA_I    ... % 17
                     OMEGA_F    ... % 18
                     LAM_1      ... % 19
                     LAM_2      ... % 20
                     LAM_3      ... % 21
                     LAM_4      ... % 22
                     LAM_5      ... % 23
                     CHI        ... % 24
                     SG         };  % 25

model.steady1 = {   kss1        ... % 1
                    Xss1        ... % 2
                    Kss1        ... % 3
                    Hss1        ... % 4
                    qss1        ... % 5
                    Yss1        ... % 6
                    yfss1       ... % 7
                    yhss1       ... % 8
                    Gss1        ... % 9
                    Zss1        ... % 10
                    Yss1        ... % 11
                    Iss1        ... % 12
                    RER         };  % 13
                
model.steady2 = {   kss2        ... % 1
                    Xss2        ... % 2
                    Kss2        ... % 3
                    Hss2        ... % 4
                    qss2        ... % 5
                    Yss2        ... % 6
                    yfss2       ... % 7
                    yhss2       ... % 8
                    Gss2        ... % 9
                    Zss2        ... % 10
                    Yss2        ... % 11
                    Iss2        ... % 12
                    RER         };  % 13

if strcmp(GLOBAL_DISPLAY, 'on')                
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp('PARETO_SS.M: STEP 3. PARAMETERS AND STEADY-STATE VALUES')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
disp(' ')
                
model.steadynames = {   'kss '        ... % 1
                        'Xss '        ... % 2
                        'Kss '        ... % 3
                        'Hss '        ... % 4
                        'qss '        ... % 5
                        'Yss '        ... % 6
                        'yfss'        ... % 7
                        'yhss'        ... % 8
                        'Gss '        ... % 9
                        'Zss '        ... % 10
                        'Yss '        ... % 11
                        'Iss '        ... % 12
                        'RER '        };  % 13
                    
model.steadyheader = {  'Var'  ...
                        'C1' ...
                        'C2'   };

s1 = [model.steady1{:}]';
s2 = [model.steady2{:}]';
steadymat = [s1, s2];

% DISPLAY TABLE OF STEADY STATE SOLUTION
disp([model.steadyheader{1},'     ', ...
                model.steadyheader{2},'     ',...
                                    model.steadyheader{3}])
disp('=========================')
for ROW = 1:length(model.steadynames)    
  disp([model.steadynames{ROW},sprintf('  %5.2f',steadymat(ROW,:))]);
end
disp('=========================')

end

disp(['K/Y ratio =',num2str(Kss1/GDP1)])

                                
