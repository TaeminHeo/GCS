clc; clear all; close all;

%define copula family
global copula_family var1_family var2_family lambda 
copula_family = 'Clayton';
var1_family = 'Gamma';
var2_family = 'Lognormal';
lambda = 60;
times = zeros(1,30);

%load data     
train = csvread(strcat('../../data/mhw/mhw_2017_train.csv'),1,0);
test = csvread(strcat('../../data/mhw/mhw_2017_test.csv'),1,0);

for i=1:30
    start = cputime;
    
    %greedy copula segmentation
    K = 1;
    seg{1} = train;
    while K > 0
        %K
        LL_all = 0;
    for j=1:K
        n = length(seg{j});
        LL_segorig = LL_v2(seg{j});
        for k=2:n-2
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
    clear seg;

    %Joint Dist. Params. Estimation
    pd1 = @(u)gamlike([u(1),u(2)],(OptimalPeriod(:,1)-u(3)));
    params1 = fminsearch(pd1,[1,1,0]);
    cdf1 = cdf(makedist(var1_family,"a",params1(1),"b",params1(2)),OptimalPeriod(:,1)-params1(3));    
    pd2 = @(u)lognlike([u(1),u(2)],(OptimalPeriod(:,2)-u(3)));
    params2 = fminsearch(pd2,[1,1,0]);
    cdf2 = cdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),OptimalPeriod(:,2)-params2(3));     
    paramhat = copulafit(copula_family,[cdf1 cdf2]);

    times(i) = cputime - start;
end

mean(times)
std(times)

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

