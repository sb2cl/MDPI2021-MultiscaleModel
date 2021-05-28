function [out] = simulateProductivity(varargin)
%% Options manager.

batch               = getFlag(varargin, 'batch');
fedbatch            = getFlag(varargin, 'fedbatch');
continous           = getFlag(varargin, 'continous');
drawFigs            = getFlag(varargin, 'drawFigs');

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

if no_combinatorial == false
    for i = 1:length(omega_A)
        for j = 1:length(kb_A)
            for k = 1:length(ku_A)
                input{i,j,k}.p.cellModel__omega_A = omega_A(i);
                input{i,j,k}.p.cellModel__kb_A = kb_A(j);
                input{i,j,k}.p.cellModel__ku_A = ku_A(k);
            end
        end
    end
else
    for i = 1:length(omega_A)
        input{i}.p.cellModel__omega_A = omega_A(i);
        input{i}.p.cellModel__kb_A = kb_A(i);
        input{i}.p.cellModel__ku_A = ku_A(i);
    end
end

%% Simulate

if batch
    % Batch
    opts.cellModel  = 'oneProtein';
    opts.massEq     = 'interpolated';
    opts.bioreactor = 'batch';
    
    % Simulate the model.
    out.batch = simulateExpCond(opts,input,'debugProductivity',false);
    
end

if fedbatch
    % Fed-Batch
    opts.cellModel  = 'oneProtein';
    opts.massEq     = 'interpolated';
    opts.bioreactor = 'fedbatch';
    
    % Simulate the model.
    out.fedbatch = simulateExpCond(opts,input);
    
end

if continous
    % Continous
    opts.cellModel  = 'oneProtein';
    opts.massEq     = 'interpolated';
    opts.bioreactor = 'continous';
    
    % Simulate the model.
    out.continous = simulateExpCond(opts,input,'debugProductivity',false);
    
end

%% Plot resutls

if batch && drawFigs
    drawProductivityPlots('batch',out.batch,varargin{:});
end

if fedbatch && drawFigs
    drawProductivityPlots('fedbatch',out.fedbatch,varargin{:});
end

if continous && drawFigs
    drawProductivityPlots('continous',out.continous,varargin{:});
end

end

%% Aux functions.

function [] = drawProductivityPlots(name, simData, varargin)

drawTrajectories = getFlag(varargin, 'drawTrajectories');
yPos             = getOption(varargin, 'yPos', 0);

if drawTrajectories
    for i = 1:n
        for j = 1:n
            simData.sp.plotAllByNamespace(batch.simOut{i,j});
        end
    end
end

figure('Position',[0 0+yPos 350 280],'Name',name);
drawFig(simData, simData.titer, 'Titer $[g \cdot L^{-1}]$', [name '_titer'], varargin{:});

figure('Position',[500 0+yPos 350 280],'Name',name);
drawFig(simData, simData.productivity, 'Productivity $[g \cdot h^{-1} \cdot L^{-1}]$', [name '_productivity'], varargin{:});

figure('Position',[1000 0+yPos 350 280],'Name',name);
drawFig(simData, simData.yield, 'Yield $[g \cdot g^{-1}]$', [name '_yield'], varargin{:});

figure('Position',[1500 0+yPos 350 280],'Name',name);
drawFig(simData, simData.muMean, '$\mu$ $[min^{-1}]$', [name '_muMean'], varargin{:});

end

function [] = drawFig(simOut, z, figTitle, filename, varargin)
%% DRAWFIG
%
% param: title
%      : mode
%      : filename
%      : drawSubplot
%      : saveFigs
%
% return:

saveFigs = getFlag(varargin, 'saveFigs');
axisXYZ  = getOption(varargin, 'axisXYZ','omega_kb');
levelZ   = getOption(varargin, 'levelZ',1);

switch(axisXYZ)
    case 'omega_kb'
        X = simOut.omega_A;
        X = squeeze(X(:,:,levelZ));
        Y = simOut.kb_A;
        Y = squeeze(Y(:,:,levelZ));
        Z = z;
        Z = squeeze(Z(:,:,levelZ));
        
    case 'omega_ku'
        X = simOut.omega_A;
        X = squeeze(X(:,levelZ,:));
        Y = simOut.ku_A;
        Y = squeeze(Y(:,levelZ,:));
        Z = z;
        Z = squeeze(Z(:,levelZ,:));
        
    case 'kb_omega'
        X = simOut.kb_A;
        X = squeeze(X(:,:,levelZ));
        Y = simOut.omega_A;
        Y = squeeze(Y(:,:,levelZ));
        Z = z;
        Z = squeeze(Z(:,:,levelZ));
        
    case 'ku_kb'
        X = simOut.ku_A;
        X = squeeze(X(levelZ,:,:));
        Y = simOut.kb_A;
        Y = squeeze(Y(levelZ,:,:));
        Z = z;
        Z = squeeze(Z(levelZ,:,:));
end

contourf(X,Y,Z, 8, '--');
colorbar;

title(figTitle,'interpreter','latex');
%xlabel('Transcription rate $N_A w_A$ $[min^{-1}]$','interpreter','latex');
xlabel('$N_A w_A$ $[min^{-1}]$','interpreter','latex');
%ylabel('Association rate RBS-ribosome $k_b^A$ $[min^{-1} \cdot molec^{-1}$]','interpreter','latex');
ylabel('$k_b^A$ $[min^{-1} \cdot molec^{-1}$]','interpreter','latex');

%   caxis([0 1.5]);
%   set(gca,'ColorScale','log');

if saveFigs
    print(['./figs/' filename],'-depsc');
end

end % drawFig
