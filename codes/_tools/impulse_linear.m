function [IRMAT] = impulse_linear(model)

IMPSELECT = model.options.IMPSELECT;
SHOCKSELECT = model.options.SHOCKSELECT;
HORIZON = model.options.HORIZON;
DO_PLOT = model.options.DO_PLOT;
MULTIPLOT = model.options.MULTIPLOT;
x_init_mat = model.options.x_init_mat;
NSHOCK = model.options.NSHOCK;  % # exogenous states

LINETYPE = model.options.LINETYPE{:};

VARNAMES = model.out.VARNAMES;
gx = model.out.gx;
hx = model.out.hx;

ny = size(gx,1); % # controls
nx = size(hx,1); % # states
nx_endog = nx - NSHOCK;  % # endogenous states
nshock_select = size(SHOCKSELECT,2);

N_show = length(IMPSELECT);

IRMAT = zeros(HORIZON+1, ny+nx, nshock_select);

   
for shock = 1:nshock_select
    x_init = x_init_mat(SHOCKSELECT(shock),:);
    IRMAT(:,:,shock) = ir(gx,hx,x_init,HORIZON+1);
    
    if DO_PLOT == 1
        if MULTIPLOT == 1
            % Plot IRFs for each shock
            figure('name',['Shock to ',...
                            VARNAMES(ny+nx_endog+SHOCKSELECT(shock),:)])
            hold on
            %suptitle2(['Shock to ',...
            %                VARNAMES(ny+nx_endog+SHOCKSELECT(shock),:)],2)
             hold off
            for i = 1:N_show
                if rem(N_show,2)==0
                    subplot(N_show/2,2,i)
                else
                    subplot((N_show+1)/2,2,i)
                end
                    hold on
                    plot(0:HORIZON,IRMAT(:,IMPSELECT(i),shock),...
                                                    LINETYPE,'LineWidth',2)
                    plot(0:HORIZON,zeros(1,HORIZON+1),'r','LineWidth',1)
                    
                    maxshock = max(abs(IRMAT(:,...
                                ny+nx_endog+SHOCKSELECT(shock),shock)));
                    maxirfcn = max(abs(IRMAT(:,IMPSELECT(i),shock)));
                    
                    if abs(maxshock) - abs(maxirfcn) < 0.005
                        plot(0:HORIZON,IRMAT(:,...
                            ny+nx_endog+SHOCKSELECT(shock),shock),'--')
                    end
                    title(VARNAMES(IMPSELECT(i),:))
                    hold off 
                    axis tight
            end
        else
                
                               
                figure('name',['Shock to ',...
                            VARNAMES(ny+nx_endog+SHOCKSELECT(shock),:)])

                plot(0:HORIZON,IRMAT(:,IMPSELECT,shock))
                legend(VARNAMES(IMPSELECT,:))
                hold on
                plot(0:HORIZON,IRMAT(:,...
                        ny+nx_endog+SHOCKSELECT(shock),shock),'--')
                hold off
                %title(VARNAMES(ny+nx_endog+SHOCKSELECT(shock),:))
                axis tight

            
        end
    end
    
    clear x_init
    
end



