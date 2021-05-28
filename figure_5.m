
out  = simulateSubstrateSweep('omega_A', [0.1], 'kb_A',[15], 'ku_A', [6],'no-combinatorial');

%%
colors;

s = out.substrateSweep.s;
productivity = out.substrateSweep.productivity;

productivity = interp1(s,productivity,0.1:0.01:3.6,'spline');
s = 0.1:0.01:3.6;

f = figure(1);
set(f,'Position',[0,0,280,240]);
clf(f);

fontSize = 8;
axisSize = 6;

% subplot(1,2,1);

hold on;

s2 = [s, fliplr(s)];
inBetween = [productivity, fliplr(ones(size(s))*productivity(end))];
h = fill(s2, inBetween,[1, 0.95, 0.85]);
% set(h,'facealpha',.5)
set(h,'edgecolor',[1 1 1])

plot(s,productivity,'Color',myColors.blue);
plot(s,ones(size(s))*productivity(end),'k:');
plot(s(end),productivity(end),'o','Color',myColors.blue);

% grid on;
set(gca,'Fontsize',axisSize);
xlabel('Substrate $[g \cdot L^{-1}]$','interpreter','latex','FontSize',fontSize);
ylabel('Productivity $[g \cdot L^{-1} \cdot h^{-1}]$','interpreter','latex','FontSize',fontSize);
ylim([0 0.005]);


% subplot(1,2,2);
% 
% s = s./max(s);
% aux = productivity;
% aux = abs(aux./aux(end)-1);
% hold on;
% 
% % s2 = [s, fliplr(s)];
% % inBetween = [productivity, fliplr(ones(size(s))*productivity(end))];
% h = fill([s, fliplr(s)], [aux,fliplr(0*ones(size(s)))],[0.7, 0.9, 1],'EdgeColor',[1 1 1]);
% % set(h,'facealpha',.5)
% % set(h,'edgecolor')
% 
% 
% plot(s,aux,'color',myColors.blue);
% 
% % grid on;
% set(gca,'Fontsize',axisSize);
% xlabel('Normalized substrate $[adim]$','interpreter','latex','FontSize',fontSize);
% ylabel('$abs(P/P_{3.6}-1)$ \quad $[adim]$','interpreter','latex','FontSize',fontSize);
% ylim([0 0.7]);


pause(1);
print('./figs/robustness_example','-depsc');
