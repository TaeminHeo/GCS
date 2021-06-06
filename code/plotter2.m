lambda = 500
result = load(strcat('../outputs/bm_result_l_',num2str(lambda),'.csv'));
result_GCS = load(strcat('../outputs/bm_result_GCS_l_',num2str(lambda),'.csv'));

diff = 100 * (result_GCS - result) ./ abs(result_GCS);
mean_diff = mean(diff);
std_diff = std(diff);

x=1:1:6;
figure()
set(gca,'FontSize',12)
patch([x fliplr(x)], [min(diff) fliplr(max(diff))],[0.8 0.8 0.8]);
hold on
line1 = plot(mean_diff,'k');
line2 = plot(min(diff),'k');
line3 = plot(max(diff),'k');
xticks([1 2 3 4 5 6]);
xlabel('m');
ylabel('\delta_L_L (%)');
saveas(gcf,strcat('../plots/BenchmarkCCAResults_',num2str(lambda),'.png'));