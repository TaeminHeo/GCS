close all;

figure(1)
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 700, 150]);
plot([1:999],[LL_all_1,0],'k')
hold on
[argvalue1, argmax1] = max(LL_all_1);
s = scatter(argmax1,argvalue1,'r');
xlabel('k');
ylabel('$\Psi_{k\backslash k+1}$','interpreter','latex');
xlim([0 1000]);
yl = ylim;
ylim([-100,yl(2)+10]);
legend(s,strcat('max (k=',num2str(argmax1),')'),'Location','best');
%print('../../plots/Psi_1','-dpng','-r300'); 

figure(2)
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 700, 150]);

subplot(1,2,1);
plot([1:argmax1-1],[LL_all_2(1,:),0],'k');
hold on
[argvalue2, argmax2] = max(LL_all_2(1,:));
s = scatter(argmax2,argvalue2,'r');
xlabel('k');
ylabel('$\Psi_{k\backslash k+1}$','interpreter','latex');
xlim([0 argmax1-1])
legend(s,strcat('max (k=',num2str(argmax2),')'),'Location','best');
yl = ylim;

subplot(1,2,2);
plot([argmax1:999],LL_all_2(2,1:(999-argmax1+1)),'k');
xlabel('k');
ylabel('$\Psi_{k\backslash k+1}$','interpreter','latex');
xlim([argmax1 1000])
ylim(yl);
%print('../../plots/Psi_2','-dpng','-r300');