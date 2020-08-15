%Prediction Performance Assessment
%test set fitting
%marginal distribution fitting
global family
exp_pd_test = fitdist(test(:,1),'Exp');
exp_cdf_test = cdf(exp_pd_test,test(:,1));
gam_pd_test = fitdist(test(:,2),'Gamma');
gam_cdf_test = cdf(gam_pd_test,test(:,2));

exp_pd_train = fitdist(train(:,1),'Exp');
exp_cdf_train = cdf(exp_pd_train,train(:,1));
gam_pd_train = fitdist(train(:,2),'Gamma');
gam_cdf_train = cdf(gam_pd_train,train(:,2));

exp_pd_optimal = fitdist(OptimalPeriod(:,1),'Exp');
exp_cdf_optimal = cdf(exp_pd_optimal,OptimalPeriod(:,1));
gam_pd_optimal = fitdist(OptimalPeriod(:,2),'Gamma');
gam_cdf_optimal = cdf(gam_pd_optimal,OptimalPeriod(:,2));

lins_d = linspace(0,20,100);
lins_s = linspace(0,15,100);
[d,s] = meshgrid(lins_d,lins_s);
[d_train,s_train] = meshgrid(cdf(exp_pd_train,lins_d),cdf(gam_pd_train,lins_s));
[d_optimal,s_optimal] = meshgrid(cdf(exp_pd_optimal,lins_d),cdf(gam_pd_optimal,lins_s));
[d_test,s_test] = meshgrid(cdf(exp_pd_test,lins_d),cdf(gam_pd_test,lins_s));
yy_train = copulacdf(family,[d_train(:),s_train(:)],paramhat_train);
yy_optimal = copulacdf(family,[d_optimal(:),s_optimal(:)],paramhat_optimal);
yy_test = copulacdf(family,[d_test(:),s_test(:)],paramhat_test);

%copula fitting
paramhat_test = copulafit(family,[exp_cdf_test gam_cdf_test]);
paramhat_train = copulafit(family,[exp_cdf_train gam_cdf_train]);
paramhat_optimal = copulafit(family,[exp_cdf_optimal gam_cdf_optimal]);

%cumulative probability
true_prob = copulacdf(family,[exp_cdf_test,gam_cdf_test],paramhat_test);
train_prob = copulacdf(family,[cdf(exp_pd_train,test(:,1)),cdf(gam_pd_train,test(:,2))],paramhat_train);
optimal_prob = copulacdf(family,[cdf(exp_pd_optimal,test(:,1)),cdf(gam_pd_optimal,test(:,2))],paramhat_optimal);

%root mean squared error (rmse)
rmse_train = sqrt(mean((true_prob - train_prob).^2));
rmse_optimal = sqrt(mean((true_prob - optimal_prob).^2));

%mean absolute percentage error (mape)
mape_train = mean(abs((true_prob - train_prob)./true_prob));
mape_optimal = mean(abs((true_prob - optimal_prob)./true_prob));

%return periods
%marginal return period
probs = 1 - 6.591549295774648 ./ [2*12 5*12 10*12 20*12 50*12 100*12];
duration_train = icdf(exp_pd_train,probs);
duration_optimal = icdf(exp_pd_optimal,probs);
duration_test = icdf(exp_pd_test,probs);

severity_train = icdf(gam_pd_train,probs);
severity_optimal = icdf(gam_pd_optimal,probs);
severity_test = icdf(gam_pd_test,probs);

%joint retrun period
close all;

rps = [2 5 10 20 50 100];

figure;
[c1,h1] = contour(d,s,6.591549295774648 / 12 ./ (1 - d_train - s_train + reshape(yy_train,100,100)),rps,'ShowText','off','color','#EDB120');
hold on
[c2,h2] = contour(d,s,6.591549295774648 / 12 ./ (1 - d_optimal - s_optimal + reshape(yy_optimal,100,100)),rps,'ShowText','off','color','#D95319');
[c3,h3] = contour(d,s,6.591549295774648 / 12 ./ (1 - d_test - s_test + reshape(yy_test,100,100)),rps,'ShowText','off','color','#0072BD');
xlabel('Duration');
ylabel('Severity');
title({'Joint Return Periods'});
legend('train set','optimal period','test set','Location','bestoutside')
%For manual labeling
%clabel(c1,h1,'manual')
%clabel(c2,h2,'manual')
%clabel(c3,h3,'manual')

%conditional retrun period
%T_D given s = 1, 5, 7

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_train,lins_d),cdf(gam_pd_train,lins_s.*0 + 1));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_train),100,100);
excd_prob = excd_prob(1,:);
d_given_s1_train = (6.591549295774648 / 12 / (1 - cdf(gam_pd_train,1))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_train,lins_d),cdf(gam_pd_train,lins_s.*0 + 5));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_train),100,100);
excd_prob = excd_prob(1,:);
d_given_s5_train = (6.591549295774648 / 12 / (1 - cdf(gam_pd_train,5))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_train,lins_d),cdf(gam_pd_train,lins_s.*0 + 7));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_train),100,100);
excd_prob = excd_prob(1,:);
d_given_s7_train = (6.591549295774648 / 12 / (1 - cdf(gam_pd_train,7))) ./ excd_prob;




[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_optimal,lins_d),cdf(gam_pd_optimal,lins_s.*0 + 1));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_optimal),100,100);
excd_prob = excd_prob(1,:);
d_given_s1_optimal = (6.591549295774648 / 12 / (1 - cdf(gam_pd_optimal,1))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_optimal,lins_d),cdf(gam_pd_optimal,lins_s.*0 + 5));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_optimal),100,100);
excd_prob = excd_prob(1,:);
d_given_s5_optimal = (6.591549295774648 / 12 / (1 - cdf(gam_pd_optimal,5))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_optimal,lins_d),cdf(gam_pd_optimal,lins_s.*0 + 7));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_optimal),100,100);
excd_prob = excd_prob(1,:);
d_given_s7_optimal = (6.591549295774648 / 12 / (1 - cdf(gam_pd_optimal,7))) ./ excd_prob;



[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_test,lins_d),cdf(gam_pd_test,lins_s.*0 + 1));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_test),100,100);
excd_prob = excd_prob(1,:);
d_given_s1_test = (6.591549295774648 / 12 / (1 - cdf(gam_pd_test,1))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_test,lins_d),cdf(gam_pd_test,lins_s.*0 + 5));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_test),100,100);
excd_prob = excd_prob(1,:);
d_given_s5_test = (6.591549295774648 / 12 / (1 - cdf(gam_pd_test,5))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_test,lins_d),cdf(gam_pd_test,lins_s.*0 + 7));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_test),100,100);
excd_prob = excd_prob(1,:);
d_given_s7_test = (6.591549295774648 / 12 / (1 - cdf(gam_pd_test,7))) ./ excd_prob;

figure;
hold on
plot(lins_d, d_given_s1_train,':','color','#EDB120')
plot(lins_d, d_given_s5_train,'--','color','#EDB120')
plot(lins_d, d_given_s7_train,'-','color','#EDB120')
plot(lins_d, d_given_s1_optimal,':','color','#D95319')
plot(lins_d, d_given_s5_optimal,'--','color','#D95319')
plot(lins_d, d_given_s7_optimal,'-','color','#D95319')
plot(lins_d, d_given_s1_test,':','color','#0072BD')
plot(lins_d, d_given_s5_test,'--','color','#0072BD')
plot(lins_d, d_given_s7_test,'-','color','#0072BD')
p1 = plot([0 0], [0 0],':k');
p2 = plot([0 0], [0 0],'--k');
p3 = plot([0 0], [0 0],'-k');
p4 = plot([0 0], [0 0],'color','#EDB120');
p5 = plot([0 0], [0 0],'color','#D95319');
p6 = plot([0 0], [0 0],'color','#0072BD');

xlabel('Duration')
ylabel('Return Period (years)')
legend([p1 p2 p3 p4 p5 p6], {'s>=1','s>=5','s>=7','train set','optimal period','test set'},'Location','best')
saveas(gcf,'../plots/d_given_s.png');

%T_S given d = 1, 5, 7

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_train,lins_d.*0 + 1),cdf(gam_pd_train,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_train),100,100);
excd_prob = excd_prob(:,1);
s_given_d1_train = (6.591549295774648 / 12 / (1 - cdf(exp_pd_train,1))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_train,lins_d.*0 + 5),cdf(gam_pd_train,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_train),100,100);
excd_prob = excd_prob(:,1);
s_given_d5_train = (6.591549295774648 / 12 / (1 - cdf(exp_pd_train,5))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_train,lins_d.*0 + 7),cdf(gam_pd_train,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_train),100,100);
excd_prob = excd_prob(:,1);
s_given_d7_train = (6.591549295774648 / 12 / (1 - cdf(exp_pd_train,7))) ./ excd_prob;




[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_optimal,lins_d.*0 + 1),cdf(gam_pd_optimal,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_optimal),100,100);
excd_prob = excd_prob(:,1);
s_given_d1_optimal = (6.591549295774648 / 12 / (1 - cdf(exp_pd_optimal,1))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_optimal,lins_d.*0 + 5),cdf(gam_pd_optimal,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_optimal),100,100);
excd_prob = excd_prob(:,1);
s_given_d5_optimal = (6.591549295774648 / 12 / (1 - cdf(exp_pd_optimal,5))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_optimal,lins_d.*0 + 7),cdf(gam_pd_optimal,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_optimal),100,100);
excd_prob = excd_prob(:,1);
s_given_d7_optimal = (6.591549295774648 / 12 / (1 - cdf(exp_pd_optimal,7))) ./ excd_prob;



[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_test,lins_d.*0 + 1),cdf(gam_pd_test,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_test),100,100);
excd_prob = excd_prob(:,1);
s_given_d1_test = (6.591549295774648 / 12 / (1 - cdf(exp_pd_test,1))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_test,lins_d.*0 + 5),cdf(gam_pd_test,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_test),100,100);
excd_prob = excd_prob(:,1);
s_given_d5_test = (6.591549295774648 / 12 / (1 - cdf(exp_pd_test,5))) ./ excd_prob;

[d_tmp,s_tmp] = meshgrid(cdf(exp_pd_test,lins_d.*0 + 7),cdf(gam_pd_test,lins_s));
excd_prob = 1 - d_tmp - s_tmp + reshape(copulacdf(family,[d_tmp(:),s_tmp(:)],paramhat_test),100,100);
excd_prob = excd_prob(:,1);
s_given_d7_test = (6.591549295774648 / 12 / (1 - cdf(exp_pd_test,7))) ./ excd_prob;

figure;
hold on
plot(lins_s, s_given_d1_train,':','color','#EDB120')
plot(lins_s, s_given_d5_train,'--','color','#EDB120')
plot(lins_s, s_given_d7_train,'-','color','#EDB120')
plot(lins_s, s_given_d1_optimal,':','color','#D95319')
plot(lins_s, s_given_d5_optimal,'--','color','#D95319')
plot(lins_s, s_given_d7_optimal,'-','color','#D95319')
plot(lins_s, s_given_d1_test,':','color','#0072BD')
plot(lins_s, s_given_d5_test,'--','color','#0072BD')
plot(lins_s, s_given_d7_test,'-','color','#0072BD')
p1 = plot([0 0], [0 0],':k');
p2 = plot([0 0], [0 0],'--k');
p3 = plot([0 0], [0 0],'-k');
p4 = plot([0 0], [0 0],'color','#EDB120');
p5 = plot([0 0], [0 0],'color','#D95319');
p6 = plot([0 0], [0 0],'color','#0072BD');

xlabel('Severity')
ylabel('Return Period (years)')
legend([p1 p2 p3 p4 p5 p6], {'d>=1','d>=5','d>=7','train set','optimal period','test set'},'Location','best')
saveas(gcf,'../plots/s_given_d.png');

