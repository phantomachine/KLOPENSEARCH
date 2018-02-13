function [meul_2,keul_2] = nb_graph_steady(xss,A,B,C,ETA,SIGMA,KAPPA,ALPHA,...
                                 BETA,GAMA,...
                                 DELTA,ZBAR,EPSILON,THETA,THETA_B,b,...
                                    OMEGA_I,OMEGA_F,PHI,TAU_X,TAU_K,TAU_H)
% NB_GRAPH_STEADY.M
% -------------------------------------------------------------------------
% Graph the steady state Money and Capital Euler equation residual
% functions, as functions of q and k, given parameterization.
%
% NOTE: You need ot have run NB_SS.M first to get some values initialized.
% -------------------------------------------------------------------------
% (c) 2010 - , Timothy Kam. Email: mortheus@gmail.com
%
% See also NB_RUN, NB_EQ, NB_SS, SSMAPSTATIC_NASHBARG, FSOLVE,
%          U_X, U_Q, G_NASH_Q, GAMMA_NASH, GAMMA_PROP
% -------------------------------------------------------------------------

qss = xss(1);
kss = xss(2);

x = [linspace(0.01,10,1000)', linspace(0.01,10,1000)'];

nx = size(x,1);
meul_1 = zeros(nx,1);
meul_2 = zeros(nx,1);
keul_1 = zeros(nx,1);
keul_2 = zeros(nx,1);

for i  = 1:size(x,1)
    [meul_1(i), keul_2(i)] = keulernb([qss, x(i,2)],A,B,C,...
                                            ETA,SIGMA,KAPPA,ALPHA,...
                                            BETA,GAMA,...
                                            DELTA,ZBAR,EPSILON,...
                                            THETA,THETA_B,b,...
                                    OMEGA_I,OMEGA_F,PHI,TAU_X,TAU_K,TAU_H);
                                
    [meul_2(i), keul_1(i)] = keulernb([x(i,1),kss],A,B,C,...
                                            ETA,SIGMA,KAPPA,ALPHA,...
                                            BETA,GAMA,...
                                            DELTA,ZBAR,EPSILON,...
                                            THETA,THETA_B,b,...
                                    OMEGA_I,OMEGA_F,PHI,TAU_X,TAU_K,TAU_H);
end

figure('name','Euler Residual Functions')

    subplot(2,1,1)
        hold on
            plot(x(:,1), zeros(length(x(:,1)),1), 'r', ...
                                    x(:,1), keul_2, 'k',...
                                                    'LineWidth',1.5);
            plot(kss,0,'o','MarkerFaceColor','r');
        hold off
        grid on
        ylabel('Capital Euler Residual, \zeta (q,k)')
        xlabel('k')

    subplot(2,1,2)
        hold on
            plot(x(:,1), zeros(length(x(:,1)),1), 'r', ...
                                    x(:,2), meul_2, 'k',...
                                                    'LineWidth',1.5);
            plot(qss,0,'o','MarkerFaceColor','r');
        hold off
        grid on
        ylabel('Money Euler Residual, \mu (q,k)')
        xlabel('q')