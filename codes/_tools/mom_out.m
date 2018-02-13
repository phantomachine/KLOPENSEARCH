function [STD_table,XCF_Y_table_sel,modelcorrcc,modelautocorr] = ...
                                mom_out(model,DATA,varshock,GLOBAL_DISPLAY)   

gx = model.out.gx;
hx = model.out.hx;

nx = size(hx,1);
ny = size(gx,1);

% Option settings
method = model.options.method;
J = model.options.J;

MOM_SELECT_x = model.options.MOM_SELECT_x;
MOM_SELECT_y = model.options.MOM_SELECT_y;

if isfield(model.options, 'CROSSCOUNTRY')
    CROSSCOUNTRY = model.options.CROSSCOUNTRY;
end

if isfield(model.options, 'XCF_INDEX')
    XCF_INDEX = model.options.XCF_INDEX;
else
    display('MOM_OUT: Setting XCF_INDEX = 1 by default!')
    XCF_INDEX = 1; % default
end

if isfield(model.options, 'Y_INDEX')
    Y_INDEX = model.options.Y_INDEX;
else
    warning('MOM_OUT:Error','You need to specify model.options.Y_INDEX')
    return
end

if ~isempty(model.options.NX1_INDEX) && ~isempty(model.options.NX2_INDEX)
    NX1_INDEX = model.options.NX1_INDEX;
    NX2_INDEX = model.options.NX2_INDEX;
end

VARNAMES = model.out.VARNAMES;

% Real data statistics: See Heathcoate and Perri, Table 2
if DATA == 1
    DATAMAT_STD = model.datastats([MOM_SELECT_y,MOM_SELECT_x], 1);
    DATAMAT_XCF_Y = model.datastats(MOM_SELECT_y, 2);
    %DATAMAT_XCF_X = model.datastats(MOM_SELECT_y, 3);
end
LAGLEAD = -J : J;

Syj = zeros(ny,ny,2*J+1);
Sxj = zeros(nx,nx,2*J+1);

% % Construct selected autocovariogram table:
ACF_table = zeros(ny+nx, 2*J+1);
XCF_table = zeros(ny, 2*J+1);
XCF_Y_table = zeros(ny, 2*J+1);
columnLabels = {1};

for j = 1 : 2*J+1
    [Syj(:,:,j), Sxj(:,:,j)] = mom(gx,hx,varshock,LAGLEAD(j),method);  

    ACF_table(:, j) = [ diag(Syj(:, :, j));
                        diag(Sxj(:, :, j)) ];

    XCF_table(:, j) = Syj(:, XCF_INDEX, j) ;
    
    XCF_Y_table(:, j) = Syj(:, Y_INDEX, j) ;
    
    columnLabels{j} = int2str(LAGLEAD(j)); 
end

% % Generate Tables:

    % Standard deviations
    STD_table = sqrt( ACF_table(:,J+1) );
    GDP_STD = STD_table(Y_INDEX);

% =======================================================================
% MOMENTS table
% =======================================================================

% All Autocorrelation functions - not used currently
for i = 1:size(ACF_table,1)
        ACF_table(i,:) = ACF_table(i,:)./ACF_table(i,J+1);
end

% Std deviations relative to GDP s.d. (%)
STD_table_sel = STD_table(MOM_SELECT_y); 

% Cross-correlation with GDP
XCF_Y_table_sel = XCF_Y_table(MOM_SELECT_y,:);

 for i = 1:size(XCF_Y_table_sel,1)
   XCF_Y_table_sel(i,:) = XCF_Y_table_sel(i,:)./(STD_table_sel(i)*GDP_STD);
 end

% Create Table: moments.tex
    NAMES_SELECT = VARNAMES(MOM_SELECT_y, :);
    var_str = cellstr(NAMES_SELECT)';
    rowLabels = cellstr(var_str);
    
if DATA == 1

    columnLabels1 = {   '(a) US Data',...
                        '(b) Model: $\frac{\% std(x)}{\% std(Y)}$',...
                        '(c) US Data', ...
                        '(d) Model: $corr(x,Y)$' };

    colmath = 0;
    rowmath = 1;

    if J == 0
    matrix2latex2([DATAMAT_STD,STD_table_sel/GDP_STD,DATAMAT_XCF_Y, XCF_Y_table_sel], ...
                    '_latex/moments.tex', ...
                    rowmath, colmath, ...
                    'rowLabels', rowLabels, ...
                    'columnLabels', columnLabels1, ...
                    'alignment', 'c', ...
                    'format', '%-6.2f',...
                    'size', 'small');
    end
    
    if strcmp(GLOBAL_DISPLAY, 'on')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('% Standard Deviation relative to GDP') 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')
    disp(['Variables', '  Model', '  US Data']);
    for n = 1:size(STD_table_sel,1)
        if ~isempty(model.options.NX1_INDEX) && ~isempty(model.options.NX2_INDEX)
            if n == NX1_INDEX || n == NX2_INDEX
            disp([NAMES_SELECT(n,:), sprintf('%-6.4f', STD_table_sel(n)), ...
                                       sprintf('%-6.4f',DATAMAT_STD(n,:))]);
            else
                disp([NAMES_SELECT(n,:), sprintf('%-6.4f', STD_table_sel(n)/GDP_STD), ...
                                       sprintf('%-6.4f',DATAMAT_STD(n,:))]);
            end
                
        else
            disp([NAMES_SELECT(n,:), sprintf('%-6.4f', STD_table_sel(n)/GDP_STD), ...
                                       sprintf('%-6.4f',DATAMAT_STD(n,:))]); 
        end
    end
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['Note: Home GDP std dev = ', sprintf('%-6.4f',GDP_STD)])
    disp(' ')
    disp(['Your chosen GDP index is model.options.Y_INDEX = ', int2str(Y_INDEX), ...
                                                          ' Check this is correct!'])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')

    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['Cross-correlation with ', VARNAMES(Y_INDEX,:)]) 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')
    disp(['Variables',  '  Model'   '  US Data']);
    disp(['Lag/Lead:  ', sprintf('  %4.0f', LAGLEAD), sprintf('  %4.0f', 0)]);
    for n = 1:size(XCF_Y_table_sel,1)       
    disp([NAMES_SELECT(n,:), sprintf('%-6.2f',XCF_Y_table_sel(n,:)),...
                                        sprintf('%-6.2f',DATAMAT_XCF_Y(n,:))]);
    end
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    end
else % DATA == 0

    columnLabels1 = {   '(a) Model: $\frac{\% std(x)}{\% std(Y)}$',...
                        '(b) Model: $corr(x,Y)$' };

    colmath = 0;
    rowmath = 1;

    if J == 0
    matrix2latex2([STD_table_sel/GDP_STD, XCF_Y_table_sel], ...
                    '_latex/moments.tex', ...
                    rowmath, colmath, ...
                    'rowLabels', rowLabels, ...
                    'columnLabels', columnLabels1, ...
                    'alignment', 'c', ...
                    'format', '%-6.4f',...
                    'size', 'small');
    end
    
    if strcmp(GLOBAL_DISPLAY, 'on')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('% Standard Deviation: relative to GDP') 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')
    disp(['Variables', '  Model']);
    for n = 1:size(STD_table_sel,1)       
    disp([NAMES_SELECT(n,:), sprintf('%-6.4f', STD_table_sel(n)/GDP_STD)]);
    end
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    disp(['Note: Home GDP std dev = ', sprintf('%-6.4f',GDP_STD)])
    disp(' ')
    disp(['Your chosen GDP index is model.options.Y_INDEX = ', int2str(Y_INDEX), ...
                                                          ' Check this is correct!'])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')

    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['Cross-correlation with ', VARNAMES(Y_INDEX,:)]) 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')
    disp(['Variables',  '  Model']);
    disp(['Lag/Lead:  ', sprintf('  %4.0f', LAGLEAD)]);
    for n = 1:size(XCF_Y_table_sel,1)       
    disp([NAMES_SELECT(n,:), sprintf('%-6.2f',XCF_Y_table_sel(n,:))]);
    end
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    end
end

%%========================================================================
% Other cross correlations
%%========================================================================

%datacorrcc = [0.58, 0.36, 0.30, 0.42]; % Data From Heathcoate-Perri Table 2

modelcorrcc = zeros(1,length(CROSSCOUNTRY));

for i = 1:length(modelcorrcc)
    modelcorrcc(i) = Syj(CROSSCOUNTRY(i,1),CROSSCOUNTRY(i,2),J+1)...
              /(STD_table(CROSSCOUNTRY(i,1))*STD_table(CROSSCOUNTRY(i,2)));
end




% if DATA == 1 && isfield(model.options, 'CROSSCOUNTRY')
%     
%     rowLabels3 = {'Data', 'Model'};
%     columnLabels3 = {'(a) $corr(Y,Y^{\ast})$', '(b) $corr(X,X^{\ast})$', ...
%                     '(c) $corr(I,I^{\ast})$', '(d) $corr(H,H^{\ast})$' };
%          
%     colmath = 0;
%     rowmath = 0;
%     
%     matrix2latex2([datacorrcc; modelcorrcc], '_latex/crosscountry.tex', ...
%                 rowmath,colmath,...
%                 'rowLabels', rowLabels3, ...
%                 'columnLabels', columnLabels3, ...
%                 'alignment', 'c', ...
%                 'format', '%-6.2f',...
%                 'size', 'small');
%             
%     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     disp('Cross-country correlations: ') 
%     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     disp(' ')
%     disp(['Correlations         ', '      Data', '      Model']);
%     for n = 1:length(modelcorrcc)       
%     disp([columnLabels3{n}, sprintf('     %-6.2f', datacorrcc(n)), ...
%                                    sprintf('     %-6.5f',modelcorrcc(n))]);
%     end
%     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     
% elseif DATA == 0 && isfield(model.options, 'CROSSCOUNTRY')
    
if isfield(model.options, 'CROSSCOUNTRY')
    
    NAMES = char(VARNAMES);
    
    rowLabels3 = {'Model'};
    
    if strcmp(GLOBAL_DISPLAY, 'on')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('Other contemporaneous correlations:')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['   (X,Z)         ', '   corr(X,Z)    ']);
    end
    
    colLabels3 = {};
    
   
    for i = 1:size(CROSSCOUNTRY,1)
        stringi = strcat('(',NAMES(CROSSCOUNTRY(i,1),:),...
                                       ',',NAMES(CROSSCOUNTRY(i,2),:),')');
        colLabels3 = [colLabels3, stringi];
        
         if strcmp(GLOBAL_DISPLAY, 'on')
            disp([stringi ,sprintf('     %-6.5f',modelcorrcc(i))]);
         end
    end
    if strcmp(GLOBAL_DISPLAY, 'on')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    end
    
    colmath = 0;
    rowmath = 0;
    
%     matrix2latex2(modelcorrcc, '_latex/crosscountry.tex', ...
%                 rowmath,colmath,...
%                 'rowLabels', rowLabels3, ...
%                 'columnLabels', colLabels3, ...
%                 'alignment', 'c', ...
%                 'format', '%-6.2f',...
%                 'size', 'small');
end


% %%========================================================================
% % Autocorrelations
% %%========================================================================
if J ~= 0
     rowLabels4 = {'Model'};
     
     VARCELL = cellstr(VARNAMES);
      
         for n = 1:length(MOM_SELECT_y)
            columnLabels4{n} = VARCELL{MOM_SELECT_y(n)};
         end
        %          
        colmath = 0;
        rowmath = 0;
        % 
        % dataautocorr = [0.85, 0.95, 0.83, 0.87, 0.91 ]; % Data From CKM Tables 5-6
        % 
        modelautocorr = ACF_table(MOM_SELECT_y,J+2)'; 
        % 
%         matrix2latex2(modelautocorr, '_latex/autocorr.tex', ...
%                         rowmath,colmath,...
%                         'rowLabels', rowLabels4, ...
%                         'columnLabels', columnLabels4, ...
%                         'alignment', 'c', ...
%                         'format', '%-6.4f',...
%                         'size', 'small');

        if strcmp(GLOBAL_DISPLAY, 'on')
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp('1st order Autocorrelations: ') 
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(' ')
        %disp(['Correlations         ', '      Model']);
        for n = 1:length(modelautocorr)       
        disp([columnLabels4{n}, sprintf('\t\t %-6.4f',modelautocorr(n))]);
        end
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        end
end % endif J~=0