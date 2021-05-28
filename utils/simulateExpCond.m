function [output] = simulateExpCond(opts, input, varargin)
%% simExpConditions
% Do a simulation starting with coherent experimental intial conditions.
%
% param: opts  Options for loading the multscale.mc model.
%        input Input parameters to simulate with.
%
% options defined by name-value pairs:
%        'tspan'     Time span for the simulation in hours.
%        'debugInit' True if we want to debug the initial conditions.
%        'debugProductivity' True if we want to debug the productivity simulations.
%
% return: output Output of the productivity simulations.

%% Options manager.

% tspan.
ind = find(ismember({varargin{1:2:end}}, 'tspan'));
if isempty(ind)
    % Default value.
    tspan = [0 +inf];
%     tspan = [0 35]*60;
else
    % User provided value.
    ind = 2 + (ind-1)*2;
    tspan = varargin{ind}*60;
end

% debugInit.
ind = find(ismember({varargin{1:2:end}}, 'debugInit'));
if isempty(ind)
    % Default value.
    debugInit = false;
else
    % User provided value.
    ind = 2 + (ind-1)*2;
    debugInit = varargin{ind};
end

% debugProductivity.
ind = find(ismember({varargin{1:2:end}}, 'debugProductivity'));
if isempty(ind)
    % Default value.
    debugProductivity = false;
else
    % User provided value.
    ind = 2 + (ind-1)*2;
    debugProductivity = varargin{ind};
end

noStopEvent = getFlag(varargin, 'noStopEvent');

%% Initialization.

% Model initialization.
fprintf('Bioreactor mode: %s\n', opts.bioreactor);
% Load the model.
m = OneModel.load('./model/multiscale.mc','opts',opts);

% Then, we initialize a SimulationClass object with the model data.
s = SimulationClass(m);

% Initialize a SimulationPlotClass object with the model data.
sp = SimulationPlotClass(m);


% Change shape of input to a 1-D vector.
[originalSize] = size(input);
input = reshape(input,[prod(originalSize),1]);
n = length(input);

%% Simulate the initial condition for the productivity experiments.

opts_init = opts;
opts_init.bioreactor = 'substrateFix';

% Load the model.
m_init = multiscale(opts_init);

% Then, we initialize a SimulationClass object with the model data.
s_init = SimulationClass(m_init);

% Initialize a SimulationPlotClass object with the model data.
sp_init = SimulationPlotClass(m_init);

for i = 1:n
    % Set the parameters of the corresponding input space.
    p = input{i}.p;
    
    % Simulate the model.
    [out_ss] = s_init.simulateSteadyState([],p,1e-6);
    
    out_ss.t = out_ss.t/60;
    
    % Save the steady state as the initial condition for the
    % productivity simulation.
    x0 = [];
    x0.cellModel__m_r = out_ss.cellModel__m_r;
    x0.cellModel__m_nr = out_ss.cellModel__m_nr;
    x0.cellModel__m_A = out_ss.cellModel__m_A;
    x0.cellModel__mu = out_ss.cellModel__mu;
    x0.cellModel__r = out_ss.cellModel__r;
    x0.bio__s = out_ss.bio__s;
    
    input{i}.x0 = x0;
    
    if i == 1
        fprintf('Initial condition progess %3.0f%%\n', i/n*100);
    else
        fprintf('\b\b\b\b\b%3.0f%%\n', i/n*100);
    end
    
    if debugInit
        figure(1);
        clf(1);
        
        % Options for the solver.
        opt = odeset('AbsTol', 1e-3, 'RelTol', 1e-6);
        
        [out] = s_init.simulate(tspan,[],p,opt);
        out.t = out.t/60;
        sp_init.plotAllStates(out);
        
        [out] = s_init.simulate(tspan,x0,p,opt);
        out.t = out.t/60;
        sp_init.plotAllStates(out);
        
        subplot(4,3,3);
        plot(out.t(end),out_ss.cellModel__mu,'o');
        
        subplot(4,3,6);
        plot(out.t(end),out_ss.cellModel__m_r,'o');
        
        subplot(4,3,9);
        plot(out.t(end),out_ss.cellModel__m_nr,'o');
        
        subplot(4,3,11);
        plot(out.t(end),out_ss.cellModel__m_A,'o');
        
        pause(0.1);
        
    end
end

%% Productivity simulation.

output = cell(n,1);
simOut = cell(n,1);

for i = 1:n
    % Set the parameters of the corresponding input space.
    try
        p = rmfield(input{i}.p,'bio__s');
    catch
        p = input{i}.p;
    end
    x0 = input{i}.x0;
    
    % Options for the solver.
    opt = odeset('AbsTol', 1e-3, 'RelTol', 1e-6);
    
    if noStopEvent == false
        if strcmp(opts.bioreactor,'batch')
            opt = odeset(opt,'Events',@(t,y,p) eventSubstrateDepletion(t,y,p,m));
        end
        
        if strcmp(opts.bioreactor,'fedbatch')
            opt = odeset(opt,'Events',@(t,y,p) eventSubstrateDepletion(t,y,p,m));
        end
        
        if strcmp(opts.bioreactor,'continous')
            opt = odeset(opt,'Events',@(t,y,p) eventFeedDepletion(t,y,p,m));
        end
    end
    
    [out] = s.simulate(tspan,x0,p,opt);
    out.t = out.t/60;
    
    % Substrate depletion time. [h]
    t_f = out.t(end);
    
    % Final amount of mass invested as protein A in the cell. [g/cell]
    mf_A = out.cellModel__m_A(end)*1e-15;
    
    % Final amount of mass of protein A removed from the bioreactor. [g]
    Mf_A = out.pro__M_A(end);
    
    % Final concentration of cells in the bioreactor. [cells/L]
    N_f = out.bio__N(end);
    
    % Initial concentration of substrate in the bioreactor. [g/L]
    S0_s = out.bio__s(1);
    
    % Final concentration of substrate in the bioreactor. [g/L]
    Sf_s = out.bio__s(end);
    
    % Substrate concentration in the feed stream. [g/L]
    Sfeed = out.bio__s_f(1);
    
    % Total mass of susbtrate removed of the bioreactor. [g]
    Sf_out = out.bio__S(end);
    
    % Intial volume of the bioreactor. [L]
    V0 = out.bio__V(1);
    
    % Final volume of the bioreactor. [L]
    Vf = out.bio__V(end);
    
    % Final volumen removed from the bioreactor. [L]
    Vf_out = out.bio__V_out(end);
    
    % Total volume feeded to the bioreactor. [L]
    Vf_feed = out.bio__V_feed(end);
    
    % Titer. [g/L]
    output{i}.titer = (mf_A*N_f*Vf + Mf_A)/(Vf + Vf_out);
    
    % Productivity. [g/h/L]
    output{i}.productivity = output{i}.titer / t_f;
    
    % Yield. [adim]
    %output{i,j}.yield = Sf_A*Vf / (S0_s*V0 - Sf_s*Vf + Sfeed*(Vf-V0));
    output{i}.yield = output{i}.titer*(Vf + Vf_out) / (S0_s*V0 - Sf_s*Vf + Sfeed*Vf_feed - Sf_out);
    
    % Mean growth rate during simulation. [min^{-1}]
    output{i}.muMean = mean(out.cellModel__mu);
    
    % Mean mu r.
    output{i}.murMean = mean(out.cellModel__mu.*out.cellModel__r);
    
    % End optical density. [OD]
    output{i}.odFinal = out.bio__N(end)*out.cellModel__m_sum(end);
    
    % Effective RBS strenght at the start of the simulation [molec^{-1}]
    output{i}.KkC0_A_start = out.cellModel__KkC0_A(1);
    
    % J_A at the start of the simulation [molec^{-1}]
    output{i}.J_A_start = out.cellModel__J_A(1);

    mid = ceil(length(out.bio__s)/2);
    
    % End biomass + substrate
    output{i}.zFinal = out.bio__s(mid) + out.bio__x(mid)*(1/out.bio__y(mid));
    
    % Initial biomass + substrate
    output{i}.zInitial = out.bio__s(1) + out.bio__x(1)*(1/out.bio__y(1));
    
    % End time. [h]
    output{i}.endTime = t_f;
    
    % Save the simulation output.
    simOut{i} = out;
    
    %fprintf('\bDone simulation %i of %i.',(i-1)*n+j,n*n);
    
    if i == 1
        fprintf('Simulation progress %3.0f%%\n', i/n*100);
    else
        fprintf('\b\b\b\b\b%3.0f%%\n', i/n*100);
    end
    
    if debugProductivity
        close all;
        sp.plotAllByNamespace(out);
        
        pause();
    end
end


%% Preparare output data for plotting.

aux = [];

for i = 1:n
    aux.omega_A(i)      = simOut{i}.cellModel__omega_A(1);
    aux.kb_A(i)         = simOut{i}.cellModel__kb_A(1);
    aux.ku_A(i)         = simOut{i}.cellModel__ku_A(1);
    aux.titer(i)        = output{i}.titer;
    aux.productivity(i) = output{i}.productivity;
    aux.yield(i)        = output{i}.yield;
    aux.muMean(i)       = output{i}.muMean;
    aux.murMean(i)      = output{i}.murMean;
    aux.KkC0_A_start(i) = output{i}.KkC0_A_start;
    aux.J_A_start(i)    = output{i}.J_A_start;
    aux.odFinal(i)      = output{i}.odFinal;
    aux.endTime(i)      = output{i}.endTime;
    aux.zFinal(i)       = output{i}.zFinal;
    aux.zInitial(i)     = output{i}.zInitial;
    aux.s(i)            = simOut{i}.bio__s(1);
end

aux.omega_A      = reshape(aux.omega_A,originalSize);
aux.kb_A         = reshape(aux.kb_A,originalSize);
aux.ku_A         = reshape(aux.ku_A,originalSize);
aux.titer        = reshape(aux.titer,originalSize);
aux.productivity = reshape(aux.productivity,originalSize);
aux.yield        = reshape(aux.yield,originalSize);
aux.muMean       = reshape(aux.muMean,originalSize);
aux.murMean      = reshape(aux.murMean,originalSize);
aux.KkC0_A_start = reshape(aux.KkC0_A_start,originalSize);
aux.J_A_start    = reshape(aux.J_A_start,originalSize);
aux.odFinal      = reshape(aux.odFinal,originalSize);
aux.endTime      = reshape(aux.endTime,originalSize);
aux.zFinal       = reshape(aux.zFinal,originalSize);
aux.zInitial     = reshape(aux.zInitial,originalSize);
aux.s            = reshape(aux.s,originalSize);

aux.simOut       = reshape(simOut,originalSize);
aux.sp           = sp;

aux.levelDiagram.inputNames{1} = 'omega_A';
aux.levelDiagram.inputNames{2} = 'kb_A';
aux.levelDiagram.inputNames{3} = 'ku_A';
aux.levelDiagram.inputNames{4} = 'KkC0_A_start';
aux.levelDiagram.inputNames{5} = 'J_A_start';

aux.levelDiagram.outputNames{1} = 'titer';
aux.levelDiagram.outputNames{2} = 'productivity';
aux.levelDiagram.outputNames{3} = 'yield';
aux.levelDiagram.outputNames{4} = 'mu mean';

aux.levelDiagram.input(:,1) = reshape(aux.omega_A,[prod(originalSize),1]);
aux.levelDiagram.input(:,2) = reshape(aux.kb_A,[prod(originalSize),1]);
aux.levelDiagram.input(:,3) = reshape(aux.ku_A,[prod(originalSize),1]);
aux.levelDiagram.input(:,4) = reshape(aux.KkC0_A_start,[prod(originalSize),1]);
aux.levelDiagram.input(:,5) = reshape(aux.J_A_start,[prod(originalSize),1]);

aux.levelDiagram.output(:,1) = reshape(aux.titer,[prod(originalSize),1]);
aux.levelDiagram.output(:,2) = reshape(aux.productivity,[prod(originalSize),1]);
aux.levelDiagram.output(:,3) = reshape(aux.yield,[prod(originalSize),1]);
aux.levelDiagram.output(:,4) = reshape(aux.muMean,[prod(originalSize),1]);

output = aux;

fprintf('Simulation finish!\n\n');

end % simulateModel
