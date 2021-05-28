%% Initialization of the model.

% [wildType oneProtein]
opts.cellModel = 'oneProtein';
% [interpolated eq_jesus fix]
opts.massEq = 'interpolated';   
% [nuFix substrateFix batch fedbatch continous]
opts.bioreactor = 'batch';

m = OneModel.load('./model/multiscale.mc','opts',opts);

%% Simulation of the model.

% Then, we initialize a SimulationClass object with the model data.
s = SimulationClass(m);

% Simulation time span.
tspan = [0 24*60];

% Parameters value we want to use in our model.
p = [];
p.cellModel__omega_A = 10;

% Initial conditions. If we dont define them, the default value is zero.
x0 = [];

% Options for the solver.
opt = odeset('AbsTol', 1e-3, 'RelTol', 1e-6);

% Simulate the model.
[out] = s.simulate(tspan,x0,p,opt);

%% Plot simulation results.

% Initialize a SimulationPlotClass object with the model data.
sp = SimulationPlotClass(m);

% Change time to hours.
out.t = out.t/60;

% close all;
sp.plotAllByNamespace(out);
