%% Figure 1 - Multiscale block diagram (show one graph of each of the scales).

close all;
clear all;

%% Simulate fedbatch mode.

% Use default values.
input{2}.p = [];
input{1}.p.cellModel__omega_A = 50;

% Batch simulation.
opts.cellModel  = 'oneProtein';
opts.massEq     = 'interpolated';
opts.bioreactor = 'batch';

out.batch = simulateExpCond(opts,input,'tspan',[0 8],'noStopEvent');

%% Plot figure 1.

close all;

fontSize = 8;
axisSize = 6;
% Bioreactor scale.
figure('Position',[0,0,235,225]);
set(gca,'Fontsize',axisSize);
hold on;
% grid on;
yyaxis right;
ylabel('Biomass $\; [g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 2]);
yticks([0:0.25:2]);
plot(out.batch.simOut{1}.t, out.batch.simOut{1}.bio__x);
plot(out.batch.simOut{2}.t, out.batch.simOut{2}.bio__x,':');
yyaxis left;
ylabel('Substrate $\; [g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 4]);
plot(out.batch.simOut{1}.t, out.batch.simOut{1}.bio__s);
plot(out.batch.simOut{2}.t, out.batch.simOut{2}.bio__s,':');
xlabel('Time $\; [h]$','interpreter','latex','FontSize',fontSize);
title('\textbf{Biorreactor model}','interpreter','latex','FontSize',fontSize);


pause(1);
print('./figs/showMultiscale_bioreactor','-depsc');

% Host scale.
figure('Position',[0,0,235,225]);
set(gca,'Fontsize',axisSize);
hold on;
% grid on;
yyaxis left;
yticks([0:0.005:0.03]);
ylabel('Growth rate $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
% set(gca,'ycolor','k') 
plot(out.batch.simOut{1}.t, out.batch.simOut{1}.cellModel__mu);
plot(out.batch.simOut{2}.t, out.batch.simOut{2}.cellModel__mu,':');

yyaxis right;
ylim([0 6000]);
ylabel('$J_{sum}(\mu,r)$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
% set(gca,'ycolor','k') 
plot(out.batch.simOut{1}.t, out.batch.simOut{1}.cellModel__J_sum);
plot(out.batch.simOut{2}.t, out.batch.simOut{2}.cellModel__J_sum,':');
xlabel('Time $\; [h]$','interpreter','latex','FontSize',fontSize);
title('\textbf{Host model}','interpreter','latex','FontSize',fontSize);
% yticks([0:750:4000]);

pause(1);
print('./figs/showMultiscale_host','-depsc');

% Synthetic circuit scale.
figure('Position',[0,0,235,225]);
set(gca,'Fontsize',axisSize);
hold on;
% grid on;
yyaxis left;
% grid on;
ylabel('Protein $A$ $\; [fg \cdot cell^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 4.5]);
plot(out.batch.simOut{1}.t, out.batch.simOut{1}.cellModel__m_A);
plot(out.batch.simOut{2}.t, out.batch.simOut{2}.cellModel__m_A,':');
yyaxis right;
% ylim([0 15]);
% yticks([0:2.5:15]);
ylabel('$N_A J_A(\mu,r)$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
plot(out.batch.simOut{1}.t, out.batch.simOut{1}.cellModel__J_A);
plot(out.batch.simOut{2}.t, out.batch.simOut{2}.cellModel__J_A,':');
xlabel('Time $\; [h]$','interpreter','latex','FontSize',fontSize);
title('\textbf{Synthetic circuit model}','interpreter','latex','FontSize',fontSize);

pause(1);
print('./figs/showMultiscale_circuit','-depsc');