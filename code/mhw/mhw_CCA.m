clc; clear all; close all;

%Get file names
itr_max = 1;

%define copula family
global copula_family var1_family var2_family lambda 
copula_family = 'Clayton';
var1_family = 'Gamma';
var2_family = 'Lognormal';
lambdas = [30,32,34,36,38,40,50,60,62,64,66,68,70,80];

for i=1:length(lambdas)
lambda = lambdas(i) %regularization parameter
do_plot = false;

%define variables
loglikelihood_GCS = zeros(itr_max,1);
loglikelihood = zeros(itr_max,1);

for itr=1:itr_max
    %load data     
    train = csvread(strcat('../../data/mhw/mhw_2017_train.csv'),1,0);
    test = csvread(strcat('../../data/mhw/mhw_2017_test.csv'),1,0);
    
    %greedy copula segmentation
    K = 1;
    seg{1} = train;
    while K > 0
        %K
        LL_all = 0;
    for j=1:K
        n = length(seg{j});
        LL_segorig = LL_v2(seg{j});
        for k=4:n-4
            LL_all(j,k) = LL_v2(seg{j}(1:k,:)) + LL_v2(seg{j}(k+1:n,:)) - LL_segorig;
        end
    end

    if max(LL_all,[],'all') > 0 
        [j_star,k_star] = find(LL_all == max(LL_all,[],'all'));
        if j_star == K
            seg{K+1} = seg{K}(k_star+1:n,:);
            seg{K} = seg{K}(1:k_star,:);
            K = K + 1;
        else
            K = -1;
        end
    else
        K = -1;
    end
    end

    %Result
    OptimalPeriod = seg{length(seg)};
    [lambda,length(OptimalPeriod)]
    OptimalPeriods{itr} = OptimalPeriod;
    Traditional{itr} = train;
    clear seg;
    
    %{
    %Joint Dist. Params. Estimation
    var1_pd_GCS = fitdist(OptimalPeriod(:,1),var1_family);
    var1_cdf_GCS = cdf(var1_pd_GCS,OptimalPeriod(:,1));
    var2_pd_GCS = fitdist(OptimalPeriod(:,2),var2_family);
    var2_cdf_GCS = cdf(var2_pd_GCS,OptimalPeriod(:,2));
    paramhat_GCS = copulafit(copula_family,[var1_cdf_GCS var2_cdf_GCS]);

    var1_pd = fitdist(train(:,1),var1_family);
    var1_cdf = cdf(var1_pd,train(:,1));
    var2_pd = fitdist(train(:,2),var2_family);
    var2_cdf = cdf(var2_pd,train(:,2));
    paramhat = copulafit(copula_family,[var1_cdf var2_cdf]);

    %Test LL
    loglikelihood_GCS(itr) = log(prod(copulapdf(copula_family,[cdf(var1_pd_GCS,test(:,1)),cdf(var2_pd_GCS,test(:,2))],paramhat_GCS),'All'));
    loglikelihood(itr) = log(prod(copulapdf(copula_family,[cdf(var1_pd,test(:,1)),cdf(var2_pd,test(:,2))],paramhat),'All'));
    %}

    %Joint Dist. Params. Estimation
    pd1_GCS = @(u)gamlike([u(1),u(2)],(OptimalPeriod(:,1)-u(3)));
    params1_GCS = fminsearch(pd1_GCS,[1,1,0]);
    cdf1_GCS = cdf(makedist(var1_family,"a",params1_GCS(1),"b",params1_GCS(2)),OptimalPeriod(:,1)-params1_GCS(3));   
    cdf1_test_GCS = cdf(makedist(var1_family,"a",params1_GCS(1),"b",params1_GCS(2)),test(:,1)-params1_GCS(3)); 
    pd2_GCS = @(u)lognlike([u(1),u(2)],(OptimalPeriod(:,2)-u(3)));
    params2_GCS = fminsearch(pd2_GCS,[1,1,0]);
    cdf2_GCS = cdf(makedist(var2_family,"mu",params2_GCS(1),"sigma",params2_GCS(2)),OptimalPeriod(:,2)-params2_GCS(3));   
    cdf2_test_GCS = cdf(makedist(var2_family,"mu",params2_GCS(1),"sigma",params2_GCS(2)),test(:,2)-params2_GCS(3));   
    paramhat_GCS = copulafit(copula_family,[cdf1_GCS cdf2_GCS]);
    
    pd1 = @(u)gamlike([u(1),u(2)],(train(:,1)-u(3)));
    params1 = fminsearch(pd1,[1,1,0]);
    cdf1 = cdf(makedist(var1_family,"a",params1(1),"b",params1(2)),train(:,1)-params1(3));   
    cdf1_test = cdf(makedist(var1_family,"a",params1(1),"b",params1(2)),test(:,1)-params1(3));   
    pd2 = @(u)lognlike([u(1),u(2)],(train(:,2)-u(3)));
    params2 = fminsearch(pd2,[1,1,0]);
    cdf2 = cdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),train(:,2)-params2(3));   
    cdf2_test = cdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),test(:,2)-params2(3));   
    paramhat = copulafit(copula_family,[cdf1 cdf2]);
    
    %Test LL
    loglikelihood_GCS(itr) = log(prod(copulapdf(copula_family,[cdf1_test_GCS,cdf2_test_GCS],paramhat_GCS),'All'));
    loglikelihood(itr) = log(prod(copulapdf(copula_family,[cdf1_test,cdf2_test],paramhat),'All'));
end

csvwrite(strcat('../../outputs/mhw/mhw_2017_result_l_',num2str(lambda),'.csv'),loglikelihood);
csvwrite(strcat('../../outputs/mhw/mhw_2017_result_GCS_l_',num2str(lambda),'.csv'),loglikelihood_GCS);
end

function loglikelihood = LL(x)
    global copula_family var1_family var2_family lambda 
    %marginal distribution fitting
    var1_pd = fitdist(x(:,1),var1_family);
    var1_cdf = cdf(var1_pd,x(:,1));
    var1_var = var(x(:,1));
    var2_pd = fitdist(x(:,2),var2_family);
    var2_cdf = cdf(var2_pd,x(:,2));
    var2_var = var(x(:,2));

    %copula fitting
    paramhat = copulafit(copula_family,[var1_cdf var2_cdf]);

    %loglikelihood
    loglikelihood = log(prod(copulapdf(copula_family,[cdf(var1_pd,x(:,1)),cdf(var2_pd,x(:,2))],paramhat),'All')) - lambda / (var1_var + var2_var);
end

function loglikelihood = LL_v2(x)
    global copula_family var1_family var2_family lambda
    MaxVal = 1000;
    options = optimset('MaxFunEvals',MaxVal,'MaxIter',MaxVal,'Display','none');
    %marginal distribution fitting
    pd1 = @(u)gamlike([u(1),u(2)],(x(:,1)-u(3)));
    params1 = fminsearch(pd1,[1,1,0],options);
    cdf1 = cdf(makedist(var1_family,"a",params1(1),"b",params1(2)),x(:,1)-params1(3));   
    var1 = var(x(:,1));

    pd2 = @(u)lognlike([u(1),u(2)],(x(:,2)-u(3)));
    params2 = fminsearch(pd2,[1,1,0],options);
    cdf2 = cdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),x(:,2)-params2(3));   
    var2 = var(x(:,2));
    
    %copula fitting
    paramhat = copulafit(copula_family,[cdf1 cdf2]);

    %loglikelihood
    loglikelihood = sum(log(copulapdf(copula_family,[cdf1,cdf2],paramhat))) - lambda / (var1 + var2);
end
