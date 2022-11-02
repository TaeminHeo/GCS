clc; clear all; close all;

lambdas = [30,32,34,36,38,40,50,60,62,64,66,68,70,80];
mean_diff = zeros(1,length(lambdas));


for i=1:length(lambdas)
    result = load(strcat('../../outputs/mhw/mhw_2017_result_l_',num2str(lambdas(i)),'.csv'));
    result_GCS = load(strcat('../../outputs/mhw/mhw_2017_result_GCS_l_',num2str(lambdas(i)),'.csv'));
    
    diff = 100 * (result_GCS - result) ./ abs(result);
    if isinf(mean(diff))
        mean_diff(i) = -100;   
    else
        mean_diff(i) = mean(diff);
    end
end

figure()
set(gca,'FontSize',12)
hold on
grid on
grid minor
plot(lambdas,mean_diff,'k');
ylim([0,90])
xlabel('\lambda')
ylabel('M (%)');
print(gcf,strcat('../../plots/mhw/mhwCCAResults.png'),'-dpng','-r300');