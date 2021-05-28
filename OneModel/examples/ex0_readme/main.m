%% 1. Model definition.

% Initialize an object of the model.
m = OneModel.load('./model/model.mc');

% Display variables and equations of the model.
disp('Variables of the model:');
m.vars
disp('Equations of the model:');
m.eqns

%% 2. Simulation.

% Initialize a SimulationClass object with the model data.
s = SimulationClass(m);

% Simulation time span.
tspan = [0 10];

% Parameters of the model.
p = []; % They are already defined in "model.mc"
% p.d1 = 10;

% Intial conditions of the model.
x0 = [];
x0.x1 = 0.000000;
x0.x2 = 0.000000;
x0.x3 = 0.000000;

% Options for the solver.
opt = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);

% Simulate the model.
[out] = s.simulate(tspan,x0,p,opt);

% Result of the simulation.
disp('Result of the simulation:');
out

%% 3. Plot simulation result.

% Initialize a SimulationPlotClass object with the model data.
sp = SimulationPlotClass(m);

% Plot the result of the simulation.
sp.plotAllStates(out);

% Save the plot.
saveas(gcf,'simulationPlot.png')

%% 4. Function that evaluates the ODEs.

% Generate an ode function of the model.
s.generateOdeFunction();
% Generate the driver script for the ode function.
s.generateDriverOdeFunction();
