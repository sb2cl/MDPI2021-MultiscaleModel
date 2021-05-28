%% Figure 2 - Productivity simulations of batch, fedbatch and continous.

close all;
clear all;

%% Simulate producitvity and add normalized values.

out = simulateProductivity('batch','fedbatch','continous');

opts.cellModel  = 'oneProtein';
opts.massEq     = 'interpolated';
opts.bioreactor = 'fedbatch';
input = {};
input{1}.p.cellModel__omega_A = 0;
input{1}.p.cellModel__kb_A = 3;
input{1}.p.cellModel__ku_A = 6;
wildType = simulateExpCond(opts,input,'tspan',[0 150]);

batch = out.batch.levelDiagram;
fedbatch = out.fedbatch.levelDiagram;
continous = out.continous.levelDiagram;

batch.outputNames = {'titer' 'titer normalized' 'productivity' 'producitivity normalized' 'yield'  'mu mean'};
batch.output = [batch.output(:,1) batch.output(:,1)./max(batch.output(:,1)) batch.output(:,2:end)];
batch.output = [batch.output(:,1:3) batch.output(:,3)./max(batch.output(:,3)) batch.output(:,4:end)];

fedbatch.outputNames = {'titer' 'titer normalized' 'productivity' 'producitivity normalized' 'yield'  'mu mean'};
fedbatch.output = [fedbatch.output(:,1) fedbatch.output(:,1)./max(fedbatch.output(:,1)) fedbatch.output(:,2:end)];
fedbatch.output = [fedbatch.output(:,1:3) fedbatch.output(:,3)./max(fedbatch.output(:,3)) fedbatch.output(:,4:end)];

continous.outputNames = {'titer' 'titer normalized' 'productivity' 'producitivity normalized' 'yield'  'mu mean'};
continous.output = [continous.output(:,1) continous.output(:,1)./max(continous.output(:,1)) continous.output(:,2:end)];
continous.output = [continous.output(:,1:3) continous.output(:,3)./max(continous.output(:,3)) continous.output(:,4:end)];

batch.outputNames = batch.outputNames([1 3 5 2 4 6]);
batch.output = batch.output(:,[1 3 5 2 4 6]);

fedbatch.outputNames = fedbatch.outputNames([1 3 5 2 4 6]);
fedbatch.output = fedbatch.output(:,[1 3 5 2 4 6]);

continous.outputNames = continous.outputNames([1 3 5 2 4 6]);
continous.output = continous.output(:,[1 3 5 2 4 6]);

%% Plot the levels diagram.

close all;

concept_batch_data = [batch.output batch.input];
concept_fedbatch_data = [fedbatch.output fedbatch.input];
concept_continous_data = [continous.output continous.input];

colors;

% Plot stylization.

fontSize = 8;
axisSize = 6;
yWildType = 1;

f = figure(1);
set(f,'Position',[0,0,560,400]);

subplot(2,3,1);
hold on;
plot(concept_batch_data(:,1),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,1),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,1),concept_continous_data(:,6),'.');
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Titer $\; [g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(2,3,4);
hold on;
plot(concept_batch_data(:,4),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,4),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,4),concept_continous_data(:,6),'.');
% plot(0,yWildType, 'ks','MarkerFaceColor',[0 0 0]);
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Normalized titer $\; [adim]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(2,3,2);
hold on;
plot(concept_batch_data(:,2),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,2),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,2),concept_continous_data(:,6),'.');
% plot(0,yWildType, 'ks','MarkerFaceColor',[0 0 0]);
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Productivity $\; [g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(2,3,5);
hold on;
plot(concept_batch_data(:,5),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,5),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,5),concept_continous_data(:,6),'.');
% plot(0,yWildType, 'ks','MarkerFaceColor',[0 0 0]);
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Normalized productivity $\; [adim]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(2,3,3);
hold on;
plot(concept_batch_data(:,3),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,3),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,3),concept_continous_data(:,6),'.');
% plot(0,yWildType, 'ks','MarkerFaceColor',[0 0 0]);
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Yield $\; [g \cdot g^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(2,3,6);
hold on;
plot(concept_batch_data(:,3)./max(concept_batch_data(:,3)),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,3)./max(concept_fedbatch_data(:,3)),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,3)./max(concept_continous_data(:,3)),concept_continous_data(:,6),'.');
% plot(0,yWildType, 'ks','MarkerFaceColor',[0 0 0]);
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Normalized Yield $\; [adim]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);


pause(1);
print('./figs/all_productivity_a','-depsc');

f = figure(2);
set(f,'Position',[0,0,370,171]);

subplot(1,2,1);
hold on;
plot(concept_batch_data(:,7),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,7),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,7),concept_continous_data(:,6),'.');
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$N_A \omega_A$ $\; [molec \cdot min^{-1} \cdot cell^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

% subplot(1,4,2);
% hold on;
% plot(concept_batch_data(:,8),concept_batch_data(:,6),'.');
% plot(concept_fedbatch_data(:,8),concept_fedbatch_data(:,6),'.');
% plot(concept_continous_data(:,8),concept_continous_data(:,6),'.');
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('$k_b^A$ $\; [min^{-1} \cdot molec^{-1}]$','interpreter','latex','FontSize',fontSize);
% % ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% 
% subplot(1,4,3);
% hold on;
% plot(concept_batch_data(:,9),concept_batch_data(:,6),'.');
% plot(concept_fedbatch_data(:,9),concept_fedbatch_data(:,6),'.');
% plot(concept_continous_data(:,9),concept_continous_data(:,6),'.');
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('$k_u^A$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% % ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);

subplot(1,2,2);
hold on;
plot(concept_batch_data(:,10),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,10),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,10),concept_continous_data(:,6),'.');
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$K^A_{C_0}(s)$ $\; [cell \cdot molec^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

pause(1);
print('./figs/all_productivity_b','-depsc');

f = figure(3);
set(f,'Position',[0,0,560,171]);

subplot(1,3,1);
hold on;
plot(concept_batch_data(:,7).*concept_batch_data(:,10),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,7).*concept_fedbatch_data(:,10),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,7).*concept_continous_data(:,10),concept_continous_data(:,6),'.');
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$N_A \omega_A K^A_{C_0}(s)$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(1,3,2);
hold on;
plot(concept_batch_data(:,11),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,11),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,11),concept_continous_data(:,6),'.');
set(gca,'Fontsize',axisSize);
grid on;
xlabel('$J_A(\mu,r,s)$ $\; [adim]$','interpreter','latex','FontSize',fontSize);
ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);

subplot(1,3,3);
hold on;
plot(concept_batch_data(:,6),concept_batch_data(:,6),'.');
plot(concept_fedbatch_data(:,6),concept_fedbatch_data(:,6),'.');
plot(concept_continous_data(:,6),concept_continous_data(:,6),'.');
plot(wildType.muMean,wildType.muMean, 'ks','MarkerFaceColor',[0 0 0]);
set(gca,'Fontsize',axisSize);
grid on;
xlabel('Mean growth rate $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.025]);
xlim([0 0.03]);


pause(1);
print('./figs/all_productivity_c','-depsc');
