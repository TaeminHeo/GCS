clc; clear all; close all;

%load data
train = csvread('../../data/Benchmark/toyproblem.csv');

%define copula family
global family lambda 
family = 'Clayton';
lambda = 100; %regularization parameter

%greedy copula segmentation
K = 1;
seg{1} = train;
while K > 0
    K
    LL_all = 0;
for j=1:K
    n = length(seg{j});
    LL_segorig = LL(seg{j});
    for k=2:n-2
        LL_all(j,k) = LL(seg{j}(1:k,:)) + LL(seg{j}(k+1:n,:)) - LL_segorig;
    end
    if K == 1
        LL_all_1 = LL_all;
    elseif K == 2
        LL_all_2 = LL_all;
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

function loglikelihood = LL(x)
    global family lambda
    %marginal distribution fitting
    pd1 = fitdist(x(:,1),'Gamma');
    cdf1 = cdf(pd1,x(:,1));
    var1 = var(x(:,1));
    pd2 = fitdist(x(:,2),'Lognormal');
    cdf2 = cdf(pd2,x(:,2));
    var2 = var(x(:,2));

    %copula fitting
    paramhat = copulafit(family,[cdf1 cdf2]);

    %loglikelihood
    loglikelihood = sum(log(copulapdf(family,[cdf1,cdf2],paramhat))) - lambda / (var1 + var2);
end

