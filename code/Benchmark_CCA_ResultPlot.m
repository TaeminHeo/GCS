close all

lambda = 10000;
result = csvread(strcat('../outputs/bm_result_l_',num2str(lambda),'.csv'));
result_GCS = csvread(strcat('../outputs/bm_result_GCS_l_',num2str(lambda),'.csv'));

diff = result_GCS(1:30,:)-result(1:30,:);
mean_diff = mean(diff);
std_diff = std(diff);

x=1:1:6;
patch([x fliplr(x)], [min(diff) fliplr(max(diff))],[0.8 0.8 0.8]);
hold on
line1 = plot(mean_diff,'k');
line2 = plot(min(diff),'k');
line3 = plot(max(diff),'k');
xticks([1 2 3 4 5 6]);
xlabel('Update Cycle');
ylabel('LL differences');
saveas(gcf,'../plots/BenchmarkCCAResults.png');
