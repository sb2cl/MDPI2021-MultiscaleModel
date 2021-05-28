function [f,c,out] = JModel(x,q)

% Mass matrix for algebraic simulations.
M = [
	0	0	0	0	
	0	0	0	0	
	0	0	1	0	
	0	0	0	1	
];

% Default initial condition value.
x0 = [
	 1.000000e-02 % cellModel__mu
	 3.500000e+02 % cellModel__r
	 1.000000e+01 % cellModel__m_r
	 1.000000e+01 % cellModel__m_nr
];

% Default parameters value.
p.cellModel__K_s = 1.802000e-01;
p.cellModel__nu_max = 1.260000e+03;
p.cellModel__m_aa = 1.826000e-07;
p.cellModel__l_e = 2.500000e+01;
p.cellModel__phi_t = 7.796000e-01;
p.cellModel__lp_r = 1.950000e+02;
p.cellModel__dm_r = 1.600000e-01;
p.cellModel__ku_r = 1.172305e+02;
p.cellModel__kb_r = 5.175400e+00;
p.cellModel__N_r = 5.500000e+01;
p.cellModel__omega_r = 5.051900e+00;
p.cellModel__ribosomeWeight = 4.500000e-03;
p.cellModel__lp_nr = 3.330000e+02;
p.cellModel__dm_nr = 2.000000e-01;
p.cellModel__ku_nr = 6.002100e+00;
p.cellModel__kb_nr = 1.116200e+01;
p.cellModel__N_nr = 1.735000e+03;
p.cellModel__omega_nr = 2.950000e-02;
p.mass__c_1 = 2.390890e+05;
p.mass__c_2 = 7.432000e+03;
p.mass__c_3 = 3.706000e+01;
p.bio__nu = 1.260000e+03;

% Use the proposed parameters of MEIGO in the model.
p.cellModel__kb_r     = x(1);
p.cellModel__ku_r     = x(2);
p.cellModel__kb_nr    = x(3);
p.cellModel__ku_nr    = x(4);
p.cellModel__omega_r  = x(5);
p.cellModel__omega_nr = x(6);
p.cellModel__phi_t    = x(7);

% Simulation time span.
tspan = [0 1000000000];

% Options for the solver.
opt = odeset('AbsTol',1e-3,'RelTol',1e-6);
opt = odeset(opt,'Mass',M);

% % Tolerance for finding the steady state.
% tol = 1e-6;
% opt = odeset(opt,'Events',@(t,y) eventSteadyState(t,y,p,tol));


for i = 1:length(q.nu)
    p.bio__nu = q.nu(i);
    [t,x] = ode15s(@(t,x) multiscaleOdeFun(t,x,p), tspan, x0, opt);
    mu(i) = x(end,1);      
end

f = sum(...
    + abs(mu' - q.mu).*[1 1 1 1 1]' ...
    );
c = [];

end