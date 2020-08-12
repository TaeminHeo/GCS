%validation
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

%copula fitting
paramhat_test = copulafit(family,[exp_cdf_test gam_cdf_test]);
paramhat_train = copulafit(family,[exp_cdf_train gam_cdf_train]);
paramhat_optimal = copulafit(family,[exp_cdf_optimal gam_cdf_optimal]);

%cumulative probability
true_prob = copulacdf(family,[exp_cdf_test,gam_cdf_test],paramhat_test);
train_prob = copulacdf(family,[cdf(exp_pd_train,test(:,1)),cdf(gam_pd_train,test(:,1))],paramhat_train);
optimal_prob = copulacdf(family,[cdf(exp_pd_optimal,test(:,1)),cdf(gam_pd_optimal,test(:,1))],paramhat_optimal);

%root mean squared error (rmse)
rmse_train = sqrt(mean((true_prob - train_prob).^2));
rmse_optimal = sqrt(mean((true_prob - optimal_prob).^2));

%mean absolute percentage error (mape)
mape_train = mean(abs((true_prob - train_prob)./true_prob));
mape_optimal = mean(abs((true_prob - optimal_prob)./true_prob));

%median absolute deviation (MAD)

