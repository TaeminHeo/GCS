close all
total = 120;
base = 20;
cycle = 10;
lambda = 150; %150

ll = load(strcat('../outputs/dr_b',num2str(base),'u',num2str(cycle),'_result_l_',num2str(lambda),'.csv'));
ll_GCS = load(strcat('../outputs/dr_b',num2str(base),'u',num2str(cycle),'_result_GCS_l_',num2str(lambda),'.csv'));

bar((0:(fix((total-base)/cycle)-1))*10,100*(ll_GCS - ll) ./ abs(ll_GCS),'Facecolor',[0.9,0.9,0.9]);
grid on
xlim([-10,100]);
xlabel('Number of additional years following the first 20 years');
ylabel('\delta_L_L (%)');
print(gcf,strcat('../plots/DroughtCCAResults_',num2str(lambda),'.png'),'-dpng','-r300');