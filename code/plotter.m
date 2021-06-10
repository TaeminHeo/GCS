close all;

figure(1)
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 700, 150]);
plot([1:999],[LL_all_1,0],'k')
hold on
[argvalue, argmax] = max(LL_all_1);
argmax
s = scatter(argmax,argvalue,'r');
xlabel('k');
ylabel('$\Psi_{k\backslash k+1}$','interpreter','latex');
ylim([0 300]);
legend(s,'max (k=722)','Location','best');
%print('../plots/Psi_1','-dpng','-r300');

figure(2)
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 700, 150]);

subplot(1,2,1);
plot([1:721],[LL_all_2(1,:),0],'k');
hold on
[argvalue, argmax] = max(LL_all_2(1,:));
argmax
s = scatter(argmax,argvalue,'r');
xlabel('k');
ylabel('$\Psi_{k\backslash k+1}$','interpreter','latex');
ylim([-20 80])
legend(s,'max (k=446)','Location','best')

subplot(1,2,2);
plot([722:999],LL_all_2(2,1:278),'k');
xlabel('k');
ylabel('$\Psi_{k\backslash k+1}$','interpreter','latex');
ylim([-20 80])
%print('../plots/Psi_2','-dpng','-r300');