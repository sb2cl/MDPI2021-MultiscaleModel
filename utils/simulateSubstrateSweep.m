function [out] = simulateSubstrateSweep(varargin)
%% Options manager.

omega_A   = getOption(varargin, 'omega_A', []);
omega_lim = getOption(varargin, 'omega_lim', [0 350]);
n_omega   = getOption(varargin, 'n_omega', 4);
kb_A      = getOption(varargin, 'kb_A', []);
n_kb      = getOption(varargin, 'n_kb', 4);
ku_A      = getOption(varargin, 'ku_A', []);
n_ku      = getOption(varargin, 'n_ku', 4);
i_max     = getOption(varargin, 'i_max', 0);

no_combinatorial = getFlag(varargin, 'no-combinatorial');


%% Definition of the input space [trascription, RBS strength].

if isempty(omega_A)
    if omega_lim(1) == 0
        % Avoid using omega_A == 0
        omega_A = arrayValuesMinMax(omega_lim(1), omega_lim(2), n_omega+1,i_max);
        omega_A = omega_A(2:end);
    else
        omega_A = arrayValuesMinMax(omega_lim(1), omega_lim(2), n_omega,i_max);
    end
end
    
if isempty(kb_A)
    kb_A    = arrayValuesMinMax(3, 15, n_kb,i_max);
end
if isempty(ku_A)
    ku_A    = arrayValuesMinMax(6, 135, n_ku,i_max);
end

s  = [0.1 0.2 0.4 0.6 0.8 1 1.5 2 2.5 3 3.6];
% s = [0.1 1 2 3.6];

input = {};

if no_combinatorial == false
    for i = 1:length(omega_A)
        for j = 1:length(kb_A)
            for k = 1:length(ku_A)
                for l = 1:length(s)
                    input{i,j,k,l}.p.cellModel__omega_A = omega_A(i);
                    input{i,j,k,l}.p.cellModel__kb_A = kb_A(j);
                    input{i,j,k,l}.p.cellModel__ku_A = ku_A(k);
                    input{i,j,k,l}.p.bio__s = s(l);
                end
            end
        end
    end
else
    for i = 1:length(omega_A)
        for j = 1:length(s)
            input{i,j}.p.cellModel__omega_A = omega_A(i);
            input{i,j}.p.cellModel__kb_A = kb_A(i);
            input{i,j}.p.cellModel__ku_A = ku_A(i);
            input{i,j}.p.bio__s = s(j);
        end
    end
end

%% Simulate

opts.cellModel  = 'oneProtein';
opts.massEq     = 'interpolated';
opts.bioreactor = 'fedbatch';

% Simulate the model.
aux = simulateExpCond(opts,input);
out.substrateSweep = aux;

% Simulate productivity.
aux = simulateProductivity('fedbatch',varargin{:});
out.productivity = aux.fedbatch;

%% Calculate the sensitivity of the substrateSweep

titerSensitivity  = [];
producSensitivity = [];
yieldSensitivity  = [];

if no_combinatorial == false
    I = size(out.substrateSweep.omega_A,1);
    J = size(out.substrateSweep.omega_A,2);
    K = size(out.substrateSweep.omega_A,3);
    
    for i = 1:I
        for j = 1:J
            for k = 1:K
                substrate  = squeeze(out.substrateSweep.s(i,j,k,:));
                titerSweep  = squeeze(out.substrateSweep.titer(i,j,k,:));
                producSweep = squeeze(out.substrateSweep.productivity(i,j,k,:));
                yieldSweep  = squeeze(out.substrateSweep.yield(i,j,k,:));
                
                % s_norm = substrate./(substrate(end));
                s_norm = substrate;
                t_norm = titerSweep./(titerSweep(end));
                p_norm = producSweep./(producSweep(end));
                y_norm = yieldSweep./(yieldSweep(end));
                
                t_robustness = trapz(s_norm,abs(t_norm-1));
                t_robustness
                
                p_robustness = trapz(s_norm,abs(p_norm-1));
                p_robustness
                
                y_robustness = trapz(s_norm,abs(y_norm-1));
                y_robustness
                
                titerSensitivity(i,j,k)  = t_robustness;
                producSensitivity(i,j,k) = p_robustness;
                yieldSensitivity(i,j,k)  = y_robustness;
                
                %             figure(1);
                %             clf(1);
                %
                %             subplot(2,3,1);
                %             plot(substrate,titerSweep);
                %             ylim([0 60]);
                %
                %             subplot(2,3,2);
                %             plot(substrate,producSweep);
                %             ylim([0 1.5]);
                %
                %             subplot(2,3,3);
                %             plot(substrate,yieldSweep);
                %             ylim([0 0.4]);
                %
                %             subplot(2,3,4);
                %             plot(substrate,abs(t_norm-1));
                %             ylim([0 1]);
                %
                %             subplot(2,3,5);
                %             plot(substrate,abs(p_norm-1));
                %             ylim([0 1]);
                %
                %             subplot(2,3,6);
                %             plot(substrate,abs(y_norm-1));
                %             ylim([0 1]);
                %
                %             pause;
                %             clc;
                
                if isnan(titerSensitivity(i,j,k))
                    titerSensitivity(i,j,k) = 0;
                end
                
                if isnan(producSensitivity(i,j,k))
                    producSensitivity(i,j,k) = 0;
                end
                
                if isnan(yieldSensitivity(i,j,k))
                    yieldSensitivity(i,j,k) = 0;
                end
            end
        end
    end
    titerSensitivity  = reshape(titerSensitivity,[I*J*K,1]);
    producSensitivity = reshape(producSensitivity,[I*J*K,1]);
    yieldSensitivity  = reshape(yieldSensitivity,[I*J*K,1]);
else
    for i = 1:size(out.substrateSweep.omega_A,1)
        
        substrate  = squeeze(out.substrateSweep.s(i,:));
        titerSweep  = squeeze(out.substrateSweep.titer(i,:));
        producSweep = squeeze(out.substrateSweep.productivity(i,:));
        yieldSweep  = squeeze(out.substrateSweep.yield(i,:));
        
        s_norm = substrate./(substrate(end));
        t_norm = titerSweep./(titerSweep(end));
        p_norm = producSweep./(producSweep(end));
        y_norm = yieldSweep./(yieldSweep(end));
        
        t_robustness = trapz(s_norm,abs(t_norm-1));
        t_robustness
        
        p_robustness = trapz(s_norm,abs(p_norm-1));
        p_robustness
        
        y_robustness = trapz(s_norm,abs(y_norm-1));
        y_robustness
        
        titerSensitivity(i)  = t_robustness;
        producSensitivity(i) = p_robustness;
        yieldSensitivity(i)  = y_robustness;
        
        if isnan(titerSensitivity(i))
            titerSensitivity(i) = 0;
        end
        
        if isnan(producSensitivity(i))
            producSensitivity(i) = 0;
        end
        
        if isnan(yieldSensitivity(i))
            yieldSensitivity(i) = 0;
        end
    end
end



%% Add the sensitivity analysis to the levelDiagram of productivity

out.productivity.levelDiagram.outputNames{7} = out.productivity.levelDiagram.outputNames{4};
out.productivity.levelDiagram.output(:,7) = out.productivity.levelDiagram.output(:,4);

out.productivity.levelDiagram.outputNames{4} = 'Titer sensitivity';
out.productivity.levelDiagram.outputNames{5} = 'Productivity sensitivity';
out.productivity.levelDiagram.outputNames{6} = 'Yield sensitivity';

out.productivity.levelDiagram.output(:,4) = titerSensitivity;
out.productivity.levelDiagram.output(:,5) = producSensitivity;
out.productivity.levelDiagram.output(:,6) = yieldSensitivity;

%% LDTool2018

% close all;
% 
% % load('./data/substrateSweep_2');
% ld = out.productivity.levelDiagram;
% 
% % input = ld.input;
% 
% % input = [input(:,1).*(input(:,2)./input(:,3))]
% 
% conceptCreate(ld.output,ld.outputNames,ld.input,ld.inputNames,'concept1');
% 
% % Bounds for Pareto front normalization
% % bounds = [upperbounds; lowerbounds]
% bounds = [
%     60 2 1 1 1 1 0.03
%     0 0 0 0 0 0 0
%     ];
% 
% % bounds = [
% %     max(ld.output)
% %     min(ld.output)
% %     ];
% % 
% % bounds(1,3:6) = [1 1 1 1];
% % bounds(2,3:6) = [0 0 0 0];
% 
% basicNorm('ld1','concept1',bounds,2);
% 
% ldDraw('ld1','concept1')
% 
% c=winter(concept1.nind);
% 
% % Ordering colors according to the values of parameter x6
% [~,idx]=sort(concept1_data(:,8));
% c2(idx,:)=c;
% ldChangeColor(ld1,c2,'concept1')

%% Plot stylization.

% close all;
% 
% c2 = [];
% c = [];

%c = copper(4);
% for i = 1:size(concept1_data,1)
%     v = concept1_data(i,8);
%     if v == 350
%         c2(i,:) = c(4,:);
%     elseif v == 262.5
%         c2(i,:) = c(3,:);
%     elseif v == 175
%         c2(i,:) = c(2,:);
%     elseif v == 87.5
%         c2(i,:) = c(1,:);
%     end
% end

% [~,idx]=sort(concept1_data(:,9));
% c=winter(concept1.nind);
% c2(idx,:)=c;
% 
% fontSize = 8;
% axisSize = 6;
% 
% f = figure(1);
% set(f,'Position',[0,0,560,400]);
% 
% subplot(2,3,1);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,1),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('Titer $\; [g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 60]);
% 
% subplot(2,3,2);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,2),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('Productivity $\; [g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 1.5]);
% 
% subplot(2,3,3);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,3),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('Yield $\; [g \cdot g^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 0.4]);
% 
% 
% subplot(2,3,4);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,4),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('Titer sensitivity $\; [1]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 1]);
% 
% subplot(2,3,5);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,5),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('Productivity sensitivity $\; [1]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 1]);
% 
% subplot(2,3,6);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,6),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('Yield sensitivity $\; [1]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 1]);
% 
% pause(1);
% print('./figs/sensitivity_high','-depsc');
% 
% f = figure(2);
% set(f,'Position',[0,0,560,171]);
% 
% subplot(1,3,1);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,8),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('$N_A \omega_A$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylabel('Growth rate $\mu$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 400]);
% 
% subplot(1,3,2);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,9),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('$k_b^A$ $\; [min^{-1} \cdot molec^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 15]);
% 
% subplot(1,3,3);
% hold on;
% for i = 1:size(concept1_data,1)
%     plot(concept1_data(i,10),concept1_data(i,7),'.','Color',c2(i,:));
% end
% set(gca,'Fontsize',axisSize);
% grid on;
% xlabel('$k_u^A$ $\; [min^{-1}]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.025]);
% xlim([0 150]);
% 
% pause(1);
% print('./figs/sensitivity_high_input','-depsc');


%% Plot brushed data

% figure(3);
% clf(3);
% 
% subplot(2,3,1);
% hold on;
% grid on;
% ylabel('Titer absolute');
% xlabel('Susbtrate');
% 
% subplot(2,3,2);
% hold on;
% grid on;
% ylabel('Productivity absolute');
% xlabel('Susbtrate');
% 
% subplot(2,3,3);
% hold on;
% grid on;
% ylabel('Yield absolute');
% xlabel('Susbtrate');
% 
% subplot(2,3,4);
% hold on;
% grid on;
% ylabel('Titer sensitivity');
% xlabel('Susbtrate');
% 
% subplot(2,3,5);
% hold on;
% grid on;
% ylabel('Productivity sensitivity');
% xlabel('Susbtrate');
% 
% subplot(2,3,6);
% hold on;
% grid on;
% ylabel('Yield sensitivity');
% xlabel('Susbtrate');
% 
% for i = 1:size(brushedData,1)
%     omega_A = brushedData(i,8);
%     kb_A    = brushedData(i,9);
%     ku_A    = brushedData(i,10);
%     
%     sim_omega_A = out.substrateSweep.omega_A;
%     sim_kb_A    = out.substrateSweep.kb_A;
%     sim_ku_A    = out.substrateSweep.ku_A;
%     
%     sim_omega_A = find(sim_omega_A == omega_A);
%     sim_kb_A    = find(sim_kb_A == kb_A);
%     sim_ku_A    = find(sim_ku_A == ku_A);
%     
%     aux   = intersect(sim_omega_A,sim_kb_A);
%     index = intersect(aux,sim_ku_A);
%     
%     susbtrate = [];
%     for j = 1:length(index)
%         susbtrate(j) = out.substrateSweep.simOut{index(j)}.bio__s(1);
%     end
%     
%     titer = out.substrateSweep.titer(index);
%     productivity = out.substrateSweep.productivity(index);
%     yield = out.substrateSweep.yield(index);
%     
%     subplot(2,3,1);
%     plot(susbtrate,titer,'Color',c2(index(1),:));
%     
%     subplot(2,3,2);
%     plot(susbtrate,productivity,'Color',c2(index(1),:));
%     
%     subplot(2,3,3);
%     plot(susbtrate,yield,'Color',c2(index(1),:));
%     
%     subplot(2,3,4);
%     plot(susbtrate,titer./max(titer),'Color',c2(index(1),:));
%     
%     subplot(2,3,5);
%     plot(susbtrate,productivity./max(productivity),'Color',c2(index(1),:));
%     
%     subplot(2,3,6);
%     plot(susbtrate,yield./max(yield),'Color',c2(index(1),:));
%     
%     % out.productivity.sp.plotAllByNamespace(out.productivity.simOut{index});
% end

end