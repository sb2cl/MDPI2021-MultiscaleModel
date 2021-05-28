function [] = simulateMassDistribution(p)
%% Default arguments.
% Model parameters to use.
if ~exist('p','var')
    p = [];
end

%% Initialization of the model.

% [wildType oneProtein]
opts.cellModel = 'wildType';
% [interpolated eq_jesus fix]
opts.massEq = 'interpolated';
% [substrateFix bioreactor]
opts.bioreactor = 'nuFix';

% m = OneModel.load('./model/multiscale.mc',opts);
m = multiscale(opts);

% Then, we initialize a SimulationClass object with the model data.
s = SimulationClass(m);

%% Simulation of the model.

% Parameters from the Bremer publication.
q = paramBremer();

n = 50;

nu = q.nu(1) + (0:n-1)*(q.nu(end) - q.nu(1))/(n-1);

OUT = [];

for i = 1:n
    % We want to use the parameters that give rise to the growth rate we
    % selected.
    p.bio__nu = nu(i);
    
    
    % Initial conditions. If we dont define them, the default value of the
    % model is used.
    x0 = [];
    
    % Options for the solver.
    opt = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
    
    % Simulate the model.
    [out] = s.simulateSteadyState(x0,p,1e-12,10,opt);
    
    if isempty(OUT)
        OUT = out;
    else
        OUT = concatStruct(OUT,out);
    end
end

%% Plot figures.

colors;

figure();

subplot(3,2,1);

hold on;
plot(q.nu,q.proteinMass,'s','Color',myColors.array{1});
plot(OUT.cellModel__nu,OUT.cellModel__m_p,'Color',myColors.array{1});
grid on;
title('Cell Mass');
legend('exp','model','Location','Best');


subplot(3,2,2);

hold on;
plot(q.nu,q.mu,'s','Color',myColors.array{1});
plot(OUT.cellModel__nu,OUT.cellModel__mu,'Color',myColors.array{1});
grid on;
title('Growth rate');
legend('exp','model','Location','Best');


subplot(3,2,3);

hold on;
plot(OUT.cellModel__nu,OUT.cellModel__m_r,'Color',myColors.array{1});
plot(OUT.cellModel__nu,OUT.cellModel__m_nr,'Color',myColors.array{2});
plot(OUT.cellModel__nu,OUT.cellModel__m_p,'Color',myColors.array{3});
plot(q.nu,q.rt.*q.r_weight,'s','Color',myColors.array{1});
plot(q.nu,q.proteinMass-q.rt.*q.r_weight,'s','Color',myColors.array{2});
plot(q.nu,q.proteinMass,'s','Color',myColors.array{3});
grid on;
title('Absolute mass distribution');
legend('model r','model nr','total mass','exp r','exp nr','exp mass','Location','Best');

subplot(3,2,4);

hold on;
plot(OUT.cellModel__nu,OUT.cellModel__m_r./OUT.cellModel__m_sum,'Color',myColors.array{1});
plot(OUT.cellModel__nu,OUT.cellModel__m_nr./OUT.cellModel__m_sum,'Color',myColors.array{2});
plot(q.nu,q.rt.*q.r_weight./q.proteinMass,'s','Color',myColors.array{1});
plot(q.nu,(q.proteinMass-q.rt.*q.r_weight)./q.proteinMass,'s','Color',myColors.array{2});
grid on;
title('Normalized mass distribution');
legend('model r','model nr','exp r','exp nr','Location','Best');

subplot(3,2,5);

hold on;
plot(q.nu,q.r,'s','Color',myColors.array{1});
plot(OUT.cellModel__nu, OUT.cellModel__r,'Color',myColors.array{1});
legend('model','exp','Location','Best');

grid on;
title('Free ribosomes');

subplot(3,2,6);

hold on;
plot(q.nu,q.rt,'s','Color',myColors.array{1});
plot(OUT.cellModel__nu, OUT.cellModel__r_t,'Color',myColors.array{1});
legend('model','exp','Location','Best');

grid on;
title('Total ribosomes');

end
