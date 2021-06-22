close all;

lambda = 10;
result = load(strcat('../outputs/bm_result_l_',num2str(lambda),'.csv'));
result_GCS = load(strcat('../outputs/bm_result_GCS_l_',num2str(lambda),'.csv'));

diff = 100 * (result_GCS - result) ./ abs(result);
mean_diff = mean(diff);
std_diff = std(diff);

x=1:1:6;
figure()
set(gca,'FontSize',12)
bar([0,100,200,300,400,500],mean_diff,'Facecolor',[0.9,0.9,0.9]);
hold on
grid on
er = errorbar([0,100,200,300,400,500],mean_diff,min(diff)-mean_diff,max(diff)-mean_diff,'k','LineStyle','none');
xlim([-80,580])
xlabel('Number of additional samples following the first 400 samples');
ylabel('\delta_L_L (%)');
legend(er,'min-max errorbars','Location','northwest');
print(gcf,strcat('../plots/BenchmarkCCAResults_',num2str(lambda),'.png'),'-dpng','-r300');