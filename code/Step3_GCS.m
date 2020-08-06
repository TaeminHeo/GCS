clc; clear all; close all;

%load data
train = csvread('../data/train.csv',1,0);
test = csvread('../data/test.csv',1,0);

%define copula family
global family 
family = 'Clayton';

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
    global family
    %marginal distribution fitting
    exp_pd = fitdist(x(:,1),'Exp');
    exp_cdf = cdf(exp_pd,x(:,1));
    gam_pd = fitdist(x(:,2),'Gamma');
    gam_cdf = cdf(gam_pd,x(:,2));

    %copula fitting
    paramhat = copulafit(family,[exp_cdf gam_cdf]);

    %loglikelihood
    loglikelihood = log(prod(copulapdf(family,[cdf(exp_pd,x(:,1)),cdf(gam_pd,x(:,2))],paramhat),'All'));
end

