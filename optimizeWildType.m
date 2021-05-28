clear all;

%% Initialization of the model.

% [wildType oneProtein]
opts.cellModel = 'wildType';
% [interpolated eq_jesus fix]
opts.massEq = 'interpolated';
% [substrateFix bioreactor]
opts.bioreactor = 'nuFix';

m = OneModel.load('./model/multiscale.mc',opts);
% m = multiscale(opts);

% Then, we initialize a SimulationClass object with the model data.
s = SimulationClass(m);

% s.generateOdeFunction();
% s.generateDriverOdeFunction();

%% MEIGO Optimization.

% Parameters from the Bremer publication.
q = paramBremer();

% We pick the values needed for the optimization. The inputs are used
% to perform the simulations, and the outputs are used to check if the
% simulation correlates with the experimental data.
aux = [];
aux.nu          = q.nu;             % Input
aux.mu          = q.mu;             % Output
aux.proteinMass = q.proteinMass;    % Output
aux.rt          = q.rt;             % Ouput
q = aux;

problem.f='JModel';
problem.x_L=[3 6 3 6 1 0.01 0];
problem.x_U=[15 135 15 135 10 0.1 1];
problem.c_L=[];
problem.c_U=[];

options.maxeval=75000;
options.local.n1=2;
options.local.n2=3;

Results=MEIGO(problem,options,'ESS',q);
delete('ess_report.mat');

vpa(Results.xbest')

p.cellModel__kb_r     = Results.xbest(1);
p.cellModel__ku_r     = Results.xbest(2);
p.cellModel__kb_nr    = Results.xbest(3);
p.cellModel__ku_nr    = Results.xbest(4);
p.cellModel__omega_r  = Results.xbest(5);
p.cellModel__omega_nr = Results.xbest(6);
p.cellModel__phi_t    = Results.xbest(7);

%% Plot result of the optimization.

simulateMassDistribution(p);
