%% Figure 3 - Sweep in kb and ku cte.

close all;
clear all;

%%

out_kb = simulateSubstrateSweep('omega_A', [1 175 350],'kb_A',[3 9 15],'ku_A',[6]);
out_ku = simulateSubstrateSweep('omega_A', [1 175 350],'kb_A',[15],'ku_A',[6 70.5 135]);

omega_A = [0.1 1 5 10 50 100 175 200 250 300 350];
out_2 = simulateSubstrateSweep('omega_A',omega_A,'n_kb',4,'n_ku',4);

wildType = simulateSubstrateSweep('omega_A', [0],'n_kb',1,'n_ku',1,'n_s');

%%

plotFigure(out_kb,out_2,wildType,'kb');

%%
plotFigure(out_ku,out_2,wildType,'ku');

%%

function [] = plotFigure(out, out_2, wildType, name)

ld = out.productivity.levelDiagram;
ld_2 = out_2.productivity.levelDiagram;

ld_color = repmat([0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.6350, 0.0780, 0.1840],3,1);
ld_size  = [9*ones(1,3) 13*ones(1,3) 18*ones(1,3)];

fontSize = 8;
axisSize = 6;
markerSize = 12;
colorAll = [0.8 0.8 0.9];

close all;

%% Plot level diagrams.

f = figure(1);
set(f,'Position',[0,0,560,400]);

subplot(2,3,1);
hold on;
plot(ld_2.output(:,1),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.output(i,1),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Titer $\; [g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
xlim([0 60]);

subplot(2,3,2);
hold on;
plot(ld_2.output(:,2),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.output(i,2),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Productivity $\; [g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
xlim([0 1.5]);

subplot(2,3,3);
hold on;
plot(ld_2.output(:,3),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.output(i,3),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Yield $\; [g \cdot g^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
xlim([0 0.4]);


subplot(2,3,4);
hold on;
plot(ld_2.output(:,4),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.output(i,4),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$RV_T$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
% xlim([0 1]);

subplot(2,3,5);
hold on;
plot(ld_2.output(:,5),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.output(i,5),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$RV_R$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
% xlim([0 1]);

subplot(2,3,6);
hold on;
plot(ld_2.output(:,6),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.output(i,6),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$RV_Y$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
% xlim([0 1]);

pause(1);
print(['./figs/sensitivity_sweep_' name],'-depsc');

%% Plot the input for the level diagrams.

f = figure(2);
set(f,'Position',[0,0,370,171]);

subplot(1,2,1);
hold on;
plot(ld_2.input(:,1),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.input(i,1),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$N_A \omega_A$ $\; [molec \cdot min^{-1} \cdot cell^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
xlim([0 400]);

subplot(1,2,2);
hold on;
plot(ld_2.input(:,4),ld_2.output(:,7),'.','Color',colorAll);
for i = size(ld.output,1):-1:1
    plot(ld.input(i,4),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$K^A_{C_0}(s_n)$ $\; [cell \cdot molec^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
xlim([0 0.3]);

% subplot(1,3,3);
% hold on;
% plot(ld_2.input(:,5),ld_2.output(:,7),'.','Color',colorAll);
% for i = size(ld.output,1):-1:1
%     plot(ld.input(i,5),ld.output(i,7),'.','MarkerSize',ld_size(i),'Color',ld_color(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('$J_A(\mu,r,s)$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 0.3]);

pause(1);
print(['./figs/sensitivity_sweep_' name '_input'],'-depsc');

%% Plot substrate sweep.

substrate = {};
titer = {};
productivity = {};
yield = {};

for i = 1:3
    for j = 1:3
        if strcmp(name,'kb')
            substrate{i,j} = squeeze(out.substrateSweep.s(i,j,1,:));
            titer{i,j} = squeeze(out.substrateSweep.titer(i,j,1,:));
            productivity{i,j} = squeeze(out.substrateSweep.productivity(i,j,1,:));
            yield{i,j} = squeeze(out.substrateSweep.yield(i,j,1,:));
        else
            substrate{i,j} = squeeze(out.substrateSweep.s(i,1,j,:));
            titer{i,j} = squeeze(out.substrateSweep.titer(i,1,j,:));
            productivity{i,j} = squeeze(out.substrateSweep.productivity(i,1,j,:));
            yield{i,j} = squeeze(out.substrateSweep.yield(i,1,j,:));
        end
    end
end

%%%

ld_color = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.6350, 0.0780, 0.1840];
% if strcmp(name,'kb')
    ld_line  = {':','--','-'};
% else
%     ld_line  = {'-','--',':'};
% end

f = figure(3);
set(f,'Position',[0,0,560,560]);


subplot(3,3,1);
hold on;
for i = 3:-1:2
    for j = 1:3
        plot((substrate{i,j}),(titer{i,j}),ld_line{j},'Color',ld_color(i,:));
    end
end
set(gca,'Fontsize',axisSize);
grid on;
ylabel('Titer $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);

subplot(3,3,2);
hold on;
for i = 2:3
    for j = 1:3
        plot((substrate{i,j}),(productivity{i,j}),ld_line{j},'Color',ld_color(i,:));
    end
end
set(gca,'Fontsize',axisSize);
grid on;
ylabel('Productivity $[g \cdot L^{-1} \cdot h^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);

subplot(3,3,3);
hold on;
for i = 2:3
    for j = 1:3
        plot((substrate{i,j}),(yield{i,j}),ld_line{j},'Color',ld_color(i,:));
    end
end
set(gca,'Fontsize',axisSize);
grid on;
ylabel('Yield $[g \cdot g^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);



subplot(3,3,4);
hold on;
for i = 1:1
    for j = 1:3
        plot((substrate{i,j}),(titer{i,j}),ld_line{j},'Color',ld_color(i,:));
    end
end
set(gca,'Fontsize',axisSize);
grid on;
ylabel('Titer $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);

subplot(3,3,5);
hold on;
for i = 1:1
    for j = 1:3
        plot(substrate{i,j},productivity{i,j},ld_line{j},'Color',ld_color(i,:));
    end
end
set(gca,'Fontsize',axisSize);
grid on;
ylabel('Productivity $[g \cdot L^{-1} \cdot h^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);

subplot(3,3,6);
hold on;
for i = 1:1
    for j = 1:3
        plot(substrate{i,j},yield{i,j},ld_line{j},'Color',ld_color(i,:));
    end
end
set(gca,'Fontsize',axisSize);
grid on;
ylabel('Yield $[g \cdot g^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);


subplot(3,3,7);
hold on; 
for i = 1:3
    for j = 1:3
        if strcmp(name,'kb')
            plot((squeeze(out.substrateSweep.s(i,j,1,:))),(squeeze(out.substrateSweep.muMean(i,j,1,:))),ld_line{j},'Color',ld_color(i,:));
        else
            plot((squeeze(out.substrateSweep.s(i,1,j,:))),(squeeze(out.substrateSweep.muMean(i,1,j,:))),ld_line{j},'Color',ld_color(i,:));
        end
    end
end
set(gca,'Fontsize',axisSize);
grid on;
plot((squeeze(wildType.substrateSweep.s)),(squeeze(wildType.substrateSweep.muMean)),'k--');
ylabel('Growth rate $[min^{-1}$]','interpreter','latex','FontSize',fontSize);
xlabel('Substrate $[g \cdot L^{-1}$]','interpreter','latex','FontSize',fontSize);

pause(1);
print(['./figs/sensitivity_sweep_' name '_substrate'],'-depsc');
end