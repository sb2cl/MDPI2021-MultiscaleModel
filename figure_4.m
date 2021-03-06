%% Figure 4 Different ways of achieving the same productivity.

% close all;
% clear all;

%% Find the interesting potins.

% We want a target productivity of 1.41 g/l/h
% With 3.6 g of substrate.
% That implies 0.000788133 1/min for growth rate
% And N_A w_A K_A = 29.168 1/molec/min

s = 3.6;

l_e = 25;
lp_A = 195;
nu_max = 1260;
K_s = 0.1802;
Emk_A = 0.62*lp_A/l_e;

nu = nu_max*s/(s+K_s);

ke = nu/l_e;

NwK = 29.168;

omega_lim = [0 350];
kb_lim = [3 15];
ku_lim = [6 135];
dm_A = 0.16;

K_min = kb_lim(1) / (ku_lim(2) + ke);
K_max = kb_lim(2) / (ku_lim(1) + ke);

out  = simulateSubstrateSweep('omega_A', [62.28 183.78 0.01], 'kb_A',[kb_lim(2) kb_lim(2) kb_lim(2)], 'ku_A', [ku_lim(1) ku_lim(1) ku_lim(1)],'no-combinatorial');

selected = {};

for i = 1:length(out.productivity.simOut)  
    selected{i}.mur = out.productivity.simOut{i}.cellModel__mu(1)*out.productivity.simOut{i}.cellModel__r(1);
    selected{i}.J_A = out.productivity.simOut{i}.cellModel__J_A(1);
end

%%

figure(1);
clf(1);

hold on;

plot(out.substrateSweep.simOut{1,1}.t,out.substrateSweep.simOut{1,1}.bio__s);
plot(out.substrateSweep.simOut{1,1}.t,out.substrateSweep.simOut{1,1}.cellModel__s);

%%

omega = [];
kb = [];
ku = [];


% low burden

kb(end+1)    = kb_lim(2);
ku(end+1)    = ku_lim(1);
omega(end+1) = fmincon(@(omega) abs(selected{1}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{1}.mur)  ),1);

kb(end+1)    = kb_lim(2);
ku(end+1)    = 20;
omega(end+1) = fmincon(@(omega) abs(selected{1}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{1}.mur)  ),1);

kb(end+1)    = kb_lim(2);
ku(end+1)    = ku_lim(2);
omega(end+1) = fmincon(@(omega) abs(selected{1}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{1}.mur)  ),1);

% high burden

kb(end+1)    = kb_lim(2);
ku(end+1)    = ku_lim(1);
omega(end+1) = fmincon(@(omega) abs(selected{2}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{2}.mur)  ),1);

kb(end+1)    = kb_lim(2);
ku(end+1)    = 20;
omega(end+1) = fmincon(@(omega) abs(selected{2}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{2}.mur)  ),1);

kb(end+1)    = kb_lim(2);
ku(end+1)    = ku_lim(2);
omega(end+1) = fmincon(@(omega) abs(selected{2}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{2}.mur)  ),1);

sep = length(omega);

% low expression

kb(end+1)    = kb_lim(2);
ku(end+1)    = ku_lim(1);
omega(end+1) = fmincon(@(omega) abs(selected{3}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{3}.mur)  ),1);

kb(end+1)    = kb_lim(2);
ku(end+1)    = 80;
omega(end+1) = fmincon(@(omega) abs(selected{3}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{3}.mur)  ),1);

kb(end+1)    = kb_lim(2);
ku(end+1)    = ku_lim(2);
omega(end+1) = fmincon(@(omega) abs(selected{3}.J_A-Emk_A*omega/(dm_A/(kb(end)/(ku(end)+ke)) + selected{3}.mur)  ),1);

%% Simulation

out  = simulateSubstrateSweep('omega_A', omega,'kb_A',kb,'ku_A',ku,'no-combinatorial');

%% Plot

colors;

fontSize = 8;
axisSize = 6;
lineWidth = 1.15;

f = figure(1);
set(f,'Position',[2000,0,380,350]);
clf(1);

lineStyle = {':',' ','--'};
lineColor = {myColors.yellow, myColors.orange};

subplot(2,2,1);
hold on;
for i = 1:sep
    plot(out.substrateSweep.s(i,:),out.substrateSweep.titer(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',lineColor{1+(i>3)});
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Titer $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 50]);

subplot(2,2,2);
hold on;
for i = 1:sep
    plot(out.substrateSweep.s(i,:),out.substrateSweep.productivity(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',lineColor{1+(i>3)});
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Productivity $[g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 1.5]);


subplot(2,2,3);
hold on;
for i = 1:sep
    plot(out.substrateSweep.s(i,:),out.substrateSweep.yield(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',lineColor{1+(i>3)});
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Yield $[g \cdot g^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.35]);

subplot(2,2,4);
hold on;
for i = 1:sep
    plot(out.substrateSweep.s(i,:),out.substrateSweep.muMean(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',lineColor{1+(i>3)});
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $[min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

pause(1);
print('./figs/sensitivity_limit_high','-depsc');


f = figure(2);
set(f,'Position',[2500,0,380,350]);
clf(2);

lineStyle = {':',' ','--'};

subplot(2,2,1);
hold on;
for i = sep+1:length(omega)
    plot(out.substrateSweep.s(i,:),out.substrateSweep.titer(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',myColors.blue);
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Titer $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.008]);

subplot(2,2,2);
hold on;
for i = sep+1:length(omega)
    plot(out.substrateSweep.s(i,:),out.substrateSweep.productivity(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',myColors.blue);
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Productivity $[g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.0005]);

subplot(2,2,3);
hold on;
for i = sep+1:length(omega)
    plot(out.substrateSweep.s(i,:),out.substrateSweep.yield(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',myColors.blue);
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Yield $[g \cdot g^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.00005]);

subplot(2,2,4);
hold on;
for i = sep+1:length(omega)
    plot(out.substrateSweep.s(i,:),out.substrateSweep.muMean(i,:),lineStyle{mod(i,3)+1},'LineWidth',1.15,'Color',myColors.blue);
end
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $[min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

pause(1);
print('./figs/sensitivity_limit_low','-depsc');



% subplot(2,4,5);
% hold on;
% for i = sep+1:length(omega)
%     plot(out.substrateSweep.s(i,:),out.substrateSweep.titer(i,:));
% end
% grid on;
% xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex');
% ylabel('Titer $[g \cdot L^{-1}]$','interpreter','latex');
% ylim([0 +inf]);
% 
% subplot(2,4,6);
% hold on;
% for i = sep+1:length(omega)
%     plot(out.substrateSweep.s(i,:),out.substrateSweep.productivity(i,:));
% end
% grid on;
% xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex');
% ylabel('Productivity $[g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex');
% ylim([0 +inf]);
% 
% 
% subplot(2,4,7);
% hold on;
% for i = sep+1:length(omega)
%     plot(out.substrateSweep.s(i,:),out.substrateSweep.yield(i,:));
% end
% grid on;
% xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex');
% ylabel('Yield $[g \cdot g^{-1}]$','interpreter','latex');
% ylim([0 +inf]);
% 
% subplot(2,4,8);
% hold on;
% for i = sep+1:length(omega)
%     t_end = [];
%     
%     for j = 1:length(out.substrateSweep.simOut)
%         t_end(j) = out.substrateSweep.simOut{i,j}.t(end);
%     end
%     plot(out.substrateSweep.s(i,:),t_end);
% end
% plot(out.substrateSweep.s(1,:),35*ones(size(out.substrateSweep.s(1,:))),'k--');
% grid on;
% xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex');
% ylabel('Completion time $[h]$','interpreter','latex');
% ylim([0 +inf]);
% ylim([0 40]);




