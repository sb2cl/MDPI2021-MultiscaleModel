%% Load the model.

clear all;

% [wildType oneProtein]
opts.cellModel = 'oneProtein';
% [interpolated eq_jesus fix]
opts.massEq = 'interpolated';   
% [substrateFix batch]
opts.bioreactor = 'batch';

m = OneModel.load('./model/multiscale.mc','opts',opts);

s = SimulationClass(m);

s.generateOdeFunction()

%% 
lc = LatexClass(m);

lc.parametersTable('./tex/parameters.tex','tab:parameters','Parameters of the model.','{\\raggedright * These parameters were re-optimized following the methods described in \\cite{nobel2020resources} to better fit the wild-type at low growth rate, since that range is the most relevant for this work. \\par}\n{\\raggedright ** Without loss of generality in the results, we choose $l_p^A$ and $d_m^A$ to be equal to the ribosomal parameters, and the range of $N_A$ and $\\omega_A$ to be in the order of ribosomal and non-ribosomal parameters. \\par}\n');
lc.variablesTable('./tex/variables.tex','tab:variables','Internal variables of the model.');
lc.statesTable('./tex/states.tex','tab:states','States and main variables of the model.');
lc.equations('./tex/equations.tex','eq:full_model');
