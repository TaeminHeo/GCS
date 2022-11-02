close all;

figure()
set(gca,'FontSize',12)
hold on
grid on
plot(lambdas,times,'k');
set(gca, 'XScale', 'log');
xlabel('\lambda')
ylabel('CPU time (s)');
%print(gcf,strcat('../../plots/DroughtTraditionalCCACPUTime.png'),'-dpng','-r300');
%print(gcf,strcat('../../plots/DroughtGCSCCACPUTime.png'),'-dpng','-r300');