% % RUN_ALL.M
% % =======================================================================
% % Script to execute sub-scripts:
% %     (1) pareto_run.m (Planner's solution - FB model)
% %     (2) compecon_run.m (Complete markets - CMFI model)
% %     (3) bondecon_run.m (Incomplete markets - CMIM model)
% % Notes:
% %     1. FB stands for "First Best" - ie. Pareto efficient allocation
% %     2. CMFI stands for "Centralized Market Full Insurance"
% %     3. CMIM stands for "Centralized Market Incomplete (Bond) Market"
% % -----------------------------------------------------------------------
% % (c) 2009 -  Timothy Kam. Email: mortheus@gmail.com
% % -----------------------------------------------------------------------

%% STEP 0: Start from clean slate
%--------------------------------------------------------------------------
clear all
close all
clc

    % PLOT MANY IN ONE FIGURE (==1), OTHERWISE (==0)
    MULTIPLOT = 0;

        fig_directory = '_figures/'; % name of figure directory
        tab_directory = '_latex/';   % LaTeX tables directory

        
%% STEP 1
%--------------------------------------------------------------------------
% Possible variables in largest set of model variables:
% WARNING: PLease check that these conincide with notation in model_eq.m
%          where 'model' refers to the economy, e.g. 'model' = 'pareto'
syms X1 X2 H1 H2 q1 q2 yh1 yf1 yf2 yh2 Y1 Y2 E I1 I2 P1 Pyh1 Pyh2 P2 ...
        Pyf2 Pyf1 RER relX ...
        K1 K2 Z1 Z2 psi1 psi2 B2
    
ymod = [ X1, X2, H1, H2, q1, q2, yh1, yf1, yf2, yh2, Y1, ...
             Y2, E, I1, I2, P1, Pyh1, Pyh2, P2, Pyf2, Pyf1, RER, relX];

xmod = [ B2, K1, K2, Z1, Z2, psi1, psi2 ];

shocksym = xmod(4:end);

shocks = {'Z','Z^{\ast}', '\psi', '\psi^{\ast}'};

varmod = [ymod, xmod];

%% STEP 2
%--------------------------------------------------------------------------
% Specify variables to plot and calculate stats:

show = [E, RER, X1, X2, q1, q2];

N_show = size(show,2);

select = zeros(1,N_show); 
for i = 1:N_show; 
    select(i) = find(varmod == show(i)); 
end

T = 30;
N_shock = 4;

% Specify collection of model name (strings):
modset = {'Pareto', 'Complete Markets', 'Incomplete Markets'};

N_model = length(modset);

%% STEP 3: Churn through all model cases
%--------------------------------------------------------------------------
    % Preallocate space for...
        % ... impulse response functions
        IRshow = zeros(T+1, N_show, N_shock, N_model);

        % ... standard deviations
        SDshow = zeros(N_show,N_model);

        % ... corr(x,GDP)
        XYshow = zeros(N_show,N_model);
        
        % ... corr(RER,E)
        REshow = zeros(1,N_model);
        
        % ... corr(RER, X/X*)
        RCshow = zeros(1,N_model);
        
        % ... corr(x, x(-1))
        ACshow = zeros(N_show,N_model);

for m = 1:N_model   
   modcase = modset{m};   
    switch modcase
       case 'Pareto' 
          disp(' ')
          disp('******************************************************** ')
          disp(['Case ', int2str(m),...
                    ': Computing and saving Pareto efficient model ...'])
          
                pareto_run;
                
                % Compute impulse response functions
                if T >  model.options.HORIZON
                    warning('MATLAB:run_all',...
                                        'T > HORIZON. Please reset T!')
                    break
                end
                
                
                if DO_IMPULSE1 == 0
                    [ IRMAT ] = impulse_linear(model);
                end
                
                falseE = ~isempty(y == 'E');
                    E_index = find(varmod == 'E');
                    if E_index == 1
                        IRMAT_1 = [zeros(T+1,falseE,2), IRMAT];
                    else
                        IRMAT_1 = [IRMAT(1:T+1,1:E_index-1,:),...
                                       zeros(T+1,falseE,size(IRMAT,3)), ...
                                               IRMAT(1:T+1,E_index:end,:)];
                    end

                
                falseP1 = ~isempty(y == 'P1');
                falsePyh1 = ~isempty(y =='Pyh1');
                falsePyh2 = ~isempty(y =='Pyh2');
                falseP2 = ~isempty(y =='P2');
                falsePyf2 = ~isempty(y =='Pyf2');
                falsePyf1 = ~isempty(y =='Pyf1');
                
                % Sub in zeros for 'P1', 'Pyh1', 'Pyh2', 'P2', 'Pyf2',
                % 'Pyf1' (Should be 6 locations)
                
                nexist = sum([falseP1,falsePyh1,falsePyh2,...
                                            falseP2,falsePyf2,falsePyf1]);
                    P1_index = find(varmod == 'P1');
                        IRMAT_1 = [ IRMAT_1(1:T+1,1:P1_index-1,:),...
                                    zeros(T+1,nexist,size(IRMAT,3)),...
                                          IRMAT_1(1:T+1,P1_index:end,:)]; 
                                      
                falseB2 = ~isempty(y == 'B2');
                    B2_index = find(varmod == 'B2');
                    if B2_index == 1
                        IRMAT_2 = [zeros(T+1,falseB2,2), IRMAT_1];
                    else
                        IRMAT_2 = [IRMAT_1(1:T+1,1:B2_index-1,:),...
                                      zeros(T+1,falseB2,size(IRMAT,3)), ...
                                            IRMAT_1(1:T+1,B2_index:end,:)];
                    end
                
                
                IRshow(:,:,1:size(IRMAT,3),m) = IRMAT_2(:,select,:);
                
                IRshow(:,:,1:size(IRMAT,3),m) = IRMAT_2(:,select,:);
                
                
                 % Compute moments   

%                     model.options.MOM_SELECT_y = 1:23;
%                     model.options.MOM_SELECT_x = [];
%                     model.options.Y_INDEX = 11;
%     
%     % Cross-country correlation selector matrix: corr(x, x*)
%     model.options.CROSSCOUNTRY = [  11, 12;                % (Y,Y*)
%                                     1,  2    ];            % (X,X*)
%                   [STD_table,XCF_Y_table,modelcorrcc,modelautocorr] ...
%                                             = mom_out(model,DATA,varshock);
                 
%                 STDshow(:,m) = STD_table();      % std deviation
%                 XCFshow(:,:,m) = XCF_Y_table();  % corr(x,GDP)
%                 CCFshow(:,:,m) = modelcorrcc;    % corr(x,z) set by CROSSCOUNTRY
%                 ACFshow(:,m) = modelautocorr();  % corr(x(t),x(t-1))
                
          
       case 'Complete Markets'
          disp(' ')
          disp('******************************************************** ')
          disp(['Case ', int2str(m),...
                    ': Computing and saving Complete Markets model ...'])
                
                compecon_run;
                
                              % Compute impulse response functions
                if T >  model.options.HORIZON
                    warning('MATLAB:run_all',...
                                        'T > HORIZON. Please reset T!')
                    break
                end
                
                if DO_IMPULSE1 == 0
                    [ IRMAT ] = impulse_linear(model);
                end
                
                falseB2 = ~isempty(y == 'B2');
                    B2_index = find(varmod == 'B2');
                    if B2_index == 1
                        IRMAT_1 = [zeros(T+1,falseB2,2), IRMAT];
                    else
                        IRMAT_1 = [IRMAT(1:T+1,1:B2_index-1,:),...
                                      zeros(T+1,falseB2,size(IRMAT,3)), ...
                                              IRMAT(1:T+1,B2_index:end,:)];
                    end
                
                
                IRshow(:,:,1:size(IRMAT,3),m) = IRMAT_1(:,select,:);
                
    
       case 'Incomplete Markets'
          disp(' ')
          disp('******************************************************** ')
          disp(['Case ', int2str(m),...
                    ': Computing and saving Incomplete Markets model ...'])
                
                bondecon_run;
                
                              % Compute impulse response functions
                if T >  model.options.HORIZON
                    warning('MATLAB:run_all',...
                                        'T > HORIZON. Please reset T!')
                    break
                end
                
                if DO_IMPULSE1 == 0
                    [ IRMAT ] = impulse_linear(model);
                end
                
                IRshow(:,:,1:size(IRMAT,3),m) = IRMAT(:,select,:);
                
      
       otherwise
          disp(' ')
          disp('******************************************************** ')
          disp('Unknown Model case.')
    end        
    close all
end

%% STEP 4: Define super-model variable listing
%--------------------------------------------------------------------------
supermodelvars;

%% STEP 5: Do graphics
%--------------------------------------------------------------------------
style = { '-.or', ...
          ':db',  ...
          '-sk',  ...
          '-pg',  ...
          '-*m'     };


if MULTIPLOT == 0
    % Construct a questdlg with 3 options
        question_str = ['This will produce N_shock x N_show =', ...
                            int2str(N_shock*N_show),' separate figures!'];
        title_str = 'run_all: ';
        [choice,killjoyflag] = choicekilljoy(question_str,title_str);

        if strcmp(choice,'Kill')
            return
        elseif strcmp(choice,'Skip')        
        elseif strcmp(choice,'Continue')
            for shoku = 1:N_shock
                for variable = 1:N_show
                    figure
                    for m = 1:N_model 
                      h = plot(0:T, IRshow(:,variable,shoku,m),...
                                                    style{m},...
                                                    'LineWidth',1.0);      
                        hold all
                    end
                    title([shocks{shoku},' \rightarrow ',...
                                            VARNAMES(select(variable),:)])
                                        
                    % Batch process and save EPS figures:
                    nowshock = char(shocksym(shoku));  % String: shock name
                    nowimpresf = char(show(variable)); % String: var name
                    
                    filename = strcat(fig_directory,...
                                        nowshock,'-to-',nowimpresf);
                    print('-depsc', filename)
                end
            end
        end
else
    
    for shoku = 1:N_shock
        figure
        for variable = 1:N_show
            if rem(N_show,2)==0
                subplot(N_show/2,2,variable)
            else
                subplot((N_show+1)/2,2,variable)
            end
            for m = 1:N_model    
               h = plot(0:T, IRshow(:,variable,shoku,m),...
                                            style{m},...
                                            'LineWidth',1.0);
                hold all
            end
            title(VARNAMES(select(variable),:))
        end
        
        % Batch process and save EPS figures:
        nowshock = char(shocksym(shoku));  % String for shock name
        
        filename = strcat(fig_directory,'irf-',nowshock);
        print('-depsc', filename)
    end
end