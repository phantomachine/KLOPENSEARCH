% DATAMASTER.M
% -------------------------------------------------------------------------
% Loads U.S. and International U.S. data (quarterly) and use HP_filter on
% some variables. Then do some summary statistics.
% -------------------------------------------------------------------------
%
% (c) 2010, T. Kam; mortheus@gmail.com
%
% See also: HPFILTER
%
% CHANGELOG: 06/09/2010, 10:00AM -- Added other OECD data, corr(RER,C/C*)
%            02/09/2010, 6:59PM

clear
clc
warning off all

PLOT_SAME = 0; % Plot in same figure
LAMBDA = 1600; % HP filter constant

EndCUT = 12;   % Remove last EndCUT # observations

% =========================================================================
%    IMPORTING OF DATA
% =========================================================================

% Data: 1975:Q1 to 2007:Q4, some missing from 2007:Q3-Q4
[data.series, text] = xlsread('alldata_fred_oecd','B2:M133');

data.series = data.series(1:end-EndCUT,:);

data.dates = text(2:end-EndCUT,1);
data.names = text(1,2:end);

% Real National accounts expenditure:
[build.series, btext] = xlsread('alldata_oecd2');

build.series = build.series(1:end-EndCUT,:);

build.dates = btext(2:end-EndCUT,1);
build.names = btext(1,2:end);

% "Rest of the World := OECD":
[row.series, btext] = xlsread('oecd_total');

row.series = row.series(1:end,:);

row.dates = btext(2:end,1);
row.names = btext(1,2:end);

% 1. CONSUMPTION
% -----------------------------------------------------------
    % FRED Consumption of Services (nominal)
    PCESV = data.series(:,1);
    
    % FRED Consumption of Non-Durables (nominal)
    PCND = data.series(:,2);
    
CONS = PCESV + PCND; % Aggregate Private Consumption (Nominal)

% 2. INVESTMENT
% -----------------------------------------------------------
    % FRED Fixed Private Investment (nominal)
    FPI = data.series(:,3);

    % FRED Consumption of Durables  (nominal)
    PCDG = data.series(:,4);
    
INVT = FPI + PCDG; % Aggregate Private Investment (Nominal)

% 3. NET EXPORTS
% -----------------------------------------------------------
    % FRED (nominal)
    NEXP = data.series(:,5);
    
% 4. DEFLATORS and REAL(OUTPUT, INVT and CONS, NEXP)
% -----------------------------------------------------------
    % FRED GDP deflator
    PDEF = data.series(:,6);
    
    % Personal Consumption Expenditures: Chain-type Price Index; Index
    % 2005=100; Q; SA; 2010-08-27
    PCECTPI = data.series(:,10);
    
    % Gross Domestic Product: Chain-type Price Index; Index 2005=100; Q;
    % SA; 2010-08-27
    GDPCTPI = data.series(:,11);
    
    % Gross Private Domestic Investment: Chain-type Price Index; Index
    % 2005=100; Q; SA; 2010-08-27
    GPDICTPI = data.series(:,12);
    
YNOM = CONS + INVT + NEXP;

%YREAL = YNOM./PDEF;
CREAL = build.series(:,1);
IREAL = build.series(:,2);
NXREAL = build.series(:,3) - build.series(:,4);

YREAL = CREAL + IREAL + NXREAL;

% Rest of the World
ROW_CREAL = row.series(:,1) - CREAL;
ROW_IREAL = row.series(:,2) - row.series(:,3) - IREAL;
ROW_NXREAL = row.series(:,4) - row.series(:,5) - NXREAL;

ROW_YREAL = ROW_CREAL + ROW_IREAL + ROW_NXREAL;

% HP_filter series (cycle) RoW's Cons, Invest., Net Exp, Output
row_CON = hpfilter(log(ROW_CREAL),LAMBDA); 
row_INV = hpfilter(log(ROW_IREAL),LAMBDA); 
row_NEX = hpfilter(ROW_NXREAL,LAMBDA); 
row_OUT = hpfilter(log(ROW_YREAL),LAMBDA); 

% 5. EMPLOYMENT DATA (see also Heathcote-Perri)
% -----------------------------------------------------------
% OECD MEI Civilian Employment Index
    AWH = data.series(:,7);

% 6. EXCHANGE RATES
% -----------------------------------------------------------
% U.S. Nominal and Real EER (IFS data)
    NEER = data.series(:,8);
    REER = data.series(:,9);
            
% =========================================================================
%    ANALYSIS OF DATA
% =========================================================================
NOBS = length(data.series);
NVAR = size(data.series,2);
N_XTICKS = 7;

% Constructed Data: From OECD Outlook Quarterly database

build.series = [YREAL, CREAL, IREAL, NXREAL, AWH, NEER, REER];
build.series(:,[1:3,5:7]) = log(build.series(:,[1:3,5:7]));

build.names = { 'Output     ',          ...
                'Consumption',          ...
                'Investment ',          ...
                'Net Exports',          ...
                'Hours      ',          ...
                'Nominal ER ',          ...
                'Real ER    '               };

NVARB = size(build.series,2);

% Do Data Plots
if PLOT_SAME == 1
    % Raw Data
    figure('Name','Raw Data Plot')
    plot(data.series)
    legend(data.names, 'Location','Best')
    xlim([0 NOBS])
    set(gca,'XTickLabel',...
                data.dates(ceil(linspace(1,NOBS,N_XTICKS))))
    
    % Constructed Data
    figure('Name','Constructed Data Plot')
    plot(build.series)
    legend(build.names, 'Location','Best')
    xlim([0 NOBS])
    set(gca,'XTickLabel',...
                data.dates(ceil(linspace(1,NOBS,N_XTICKS))))
else
    
    % Raw Data
    figure('Name','Raw Data Plot')
    for i = 1 : NVAR
        subplot((NVAR+rem(NVAR,2))/2,2,i)
        plot(data.series(:,i))
        ylabel(data.names{i})
        xlim([0 NOBS])
        set(gca,'XTickLabel',...
                data.dates(ceil(linspace(1,NOBS,N_XTICKS))))
    end
    
     % Constructed Data
    figure('Name','Constructed Data Plot')
    for i = 1 : NVARB
        subplot((NVARB+rem(NVARB,2))/2,2,i)
        plot(build.series(:,i))
        ylabel(build.names{i})
        xlim([0 NOBS])
        set(gca,'XTickLabel',...
                data.dates(ceil(linspace(1,NOBS,N_XTICKS))))
    end
end    

% =========================================================================
% HP-filter the data
% =========================================================================
HP_SELECT = 1:NVARB;
nLags = 1;
nSTDs = 2;

N_HP = length(HP_SELECT);

build.hpf = zeros(NOBS,N_HP);
build.hpf_trend = build.hpf;
build.hpf_cycle = build.hpf;

STD = zeros(N_HP,1);            %STD_bounds = zeros(N_HP,2);
ACF = zeros(N_HP,nLags+1);
XCF = zeros(N_HP,2*nLags+1);

for i = 1 : N_HP
    y_temp = build.series(:,HP_SELECT(i));
    
    % HP_filter data
    [cycle_temp, build.hpf_trend(:,i)] = hpfilter(y_temp,LAMBDA); 
    
    build.hpf_cycle(:,i) = cycle_temp;
    
    % Calculate moments (std)
    STD(i) = std(cycle_temp);
    ACF(i,:) = autocorr(cycle_temp, nLags)';
end

clear y_temp

% Do Data Plots
if PLOT_SAME == 1
        
    % Constructed Data
    figure('Name','HP-filtered Data')
    plot(build.hpf_cycle(:,:))
    legend(build.names, 'Location','Best')
    xlim([0 NOBS])
    set(gca,'XTickLabel',...
                data.dates(ceil(linspace(1,NOBS,N_XTICKS))))
else
    
    % Constructed Data
    figure('Name','HP-filtered Data')
    for i = 1 : NVARB
        subplot((NVARB+rem(NVARB,2))/2,2,i)
        plot(build.hpf_cycle(:,i),'-+r')
        ylabel(build.names{i})
        xlim([0 NOBS])
        set(gca,'XTickLabel',...
                data.dates(ceil(linspace(1,NOBS,N_XTICKS))))
    end
end    


% Build tables

Y_INDEX = 1;
CON_INDEX = 2;
INV_INDEX = 3;
NEX_INDEX = 4;
NER_INDEX = 6;
RER_INDEX = 7;

                           
                            
% Standard deviations
disp(' ')
disp('-------------------------------------------------------------------')
disp(['Variable   ', '      Std dev (Percent)*     ', ...
            ' Std relative to Output'  ])
disp('-------------------------------------------------------------------')
for i = 1 : N_HP
    if i == NEX_INDEX
        disp([build.names{i},  sprintf('\t\t %-6.2f', STD(i)*100), ...
                        sprintf('\t n.a. %-6.2f', [])])
    else
        disp([build.names{i},  sprintf('\t\t %-6.2f', STD(i)*100), ...
                        sprintf('\t\t %-6.2f', STD(i)/STD(Y_INDEX))])
    end
end
disp('-------------------------------------------------------------------')
disp('* Note: With the exception of Net Exports (in levels).')
disp('-------------------------------------------------------------------')

% 1st order autocorrelation
disp(' ')
disp('-------------------------------------------------------------------')
disp(['Variable   ', '      Autocorrelation     '])
disp('-------------------------------------------------------------------')
for i = 1 : N_HP
    
        disp([build.names{i},  sprintf('\t\t %-6.2f', ACF(i,2))])
    
end
disp('-------------------------------------------------------------------')


% Other correlations

CORR_RER_NER = crosscorr(build.hpf_cycle(:,NER_INDEX), ...
                                build.hpf_cycle(:,RER_INDEX), nLags)';

CORR_RER_NEX = crosscorr(build.hpf_cycle(:,NEX_INDEX), ...
                                build.hpf_cycle(:,RER_INDEX), nLags)'; 
                            
CORR_CON_ROWCON = crosscorr(build.hpf_cycle(:,CON_INDEX), ...
                                                    row_CON, nLags)';
                                                
CORR_Y_ROWY = crosscorr(build.hpf_cycle(:,Y_INDEX), ...
                                                    row_OUT, nLags)';
                                                
CORR_I_ROWI = crosscorr(build.hpf_cycle(:,INV_INDEX), ...
                                                    row_INV, nLags)';
CORR_NEX_ROWNEX = crosscorr(build.hpf_cycle(:,NEX_INDEX), ...
                                                    row_NEX, nLags)';

REL_CON = build.hpf_cycle(:,2) - row_CON;    % log(C/C*)

CORR_RELCON_RER = crosscorr(REL_CON, ...
                                build.hpf_cycle(:,RER_INDEX), nLags)';

disp(' ')
disp('-------------------------------------------------------------------')
disp(['Variable   ', '            Cross-country Correlation     '])
disp('-------------------------------------------------------------------')

 disp(['corr(','NER',', ','RER',')',...
                    sprintf('\t\t %-6.2f', CORR_RER_NER(:,nLags+1))])
                
 disp(['corr(','NEX',', ','RER',')',...
                    sprintf('\t\t %-6.2f', CORR_RER_NEX(:,nLags+1))])
 disp(['corr(','CON',', ','CON*',')',...
                    sprintf('\t\t %-6.2f', CORR_CON_ROWCON(:,nLags+1))])  
 disp(['corr(',' Y ',', ', ' Y* ',')',...
                    sprintf('\t\t %-6.2f', CORR_Y_ROWY(:,nLags+1))])  
 disp(['corr(',' I ',', ', ' I* ',')',...
                    sprintf('\t\t %-6.2f', CORR_I_ROWI(:,nLags+1))]) 
 disp(['corr(','NEX',', ', 'NEX*',')',...
                    sprintf('\t\t %-6.2f', CORR_NEX_ROWNEX(:,nLags+1))]) 
 disp(['corr(','RER',', ', 'C/C*',')',...
                    sprintf('\t\t %-6.2f', CORR_RELCON_RER(:,nLags+1))])               

disp('-------------------------------------------------------------------')
disp('* Note: R.o.W. is OECD aggregate less USA')
disp('-------------------------------------------------------------------')

disp(' ')
disp('-------------------------------------------------------------------')
disp(['Your data sample is from ',data.dates{1}, ' to ', data.dates{end}])
disp('-------------------------------------------------------------------')


    