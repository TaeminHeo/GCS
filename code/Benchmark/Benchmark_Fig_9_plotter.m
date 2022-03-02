clc; clear all; close all;

lambda = [1,5,10,50,100,500,1000,5000,10000];
mean_diff = zeros(1,length(lambda));
min_diff = zeros(1,length(lambda));
max_diff = zeros(1,length(lambda));

for i=1:length(lambda)
    result = load(strcat('../../outputs/bm_result_l_',num2str(lambda(i)),'.csv'));
    result_GCS = load(strcat('../../outputs/bm_result_GCS_l_',num2str(lambda(i)),'.csv'));
    diff = 100 * (result_GCS - result) ./ abs(result);
    mean_diff(i) = mean(diff,'all');
    min_diff(i) = min(mean(diff,2));
    max_diff(i) = max(mean(diff,2));
end

figure()
set(gca,'FontSize',12)
hold on
grid on
plot(lambda,mean_diff,'k','DisplayName','mean of M');
errorbar(lambda,mean_diff,min_diff-mean_diff,max_diff-mean_diff,'k','LineStyle','none','DisplayName','min-max errorbars');
set(gca, 'XScale', 'log');
xlim([0.8,12000])
xlabel('\lambda')
ylabel('M (%)');
legend();
%print(gcf,strcat('../../plots/BenchmarkCCAResults.png'),'-dpng','-r300');