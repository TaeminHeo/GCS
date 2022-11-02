clc; clear all; close all;

base = 20;
cycle = 10;
lambda = [0.001,0.01,0.1];
mean_diff = zeros(1,length(lambda));


for i=1:length(lambda)
    result = load(strcat('../../outputs/SalinitySST/ss_b',num2str(base),'u',num2str(cycle),'_result_l_',num2str(lambda(i)),'.csv'));
    result_GCS = load(strcat('../../outputs/SalinitySST/ss_b',num2str(base),'u',num2str(cycle),'_result_GCS_l_',num2str(lambda(i)),'.csv'));
    
    diff = 100 * (result_GCS - result) ./ abs(result);
    mean_diff(i) = mean(diff);    
end

%x=1:1:6;
figure()
set(gca,'FontSize',12)
hold on
grid on
grid minor
plot(lambda,mean_diff,'k');
%ylim([-10,10])
xlabel('\lambda')
ylabel('M (%)');
%print(gcf,strcat('../../plots/DroughtCCAResults.png'),'-dpng','-r300');