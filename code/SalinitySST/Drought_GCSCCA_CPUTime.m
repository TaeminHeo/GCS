clc; clear all; close all;

%Get file names
total = 120;
base = 20;
cycle = 10;
fnames = dir(strcat('../../data/b',num2str(base),'u',num2str(cycle),'_*.csv'));
fnames = {fnames.name};

%define copula family
global family lambda 
family = 'Gumbel';
lambdas = [50,100,110,120,130,140,150,160,170,180,190,200,250];
times = zeros(1,length(lambdas));

for i=1:length(lambdas)
lambda = lambdas(i); %regularization parameter
do_plot = false;

%define variables
loglikelihood_GCS = zeros(fix((total-base)/cycle),1);
start = cputime;
for itr=1:fix((total-base)/cycle)
    %load data     
    train = csvread(strcat('../../data/b',num2str(base),'u',num2str(cycle),'_',num2str(itr-1),'.csv'),1,0);
    test = csvread(strcat('../../data/b',num2str(base),'u',num2str(cycle),'_',num2str(itr-1),'_test.csv'),1,0);
    
    %greedy copula segmentation
    K = 1;
    seg{1} = train;
    while K > 0
        %K
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
    OptimalPeriods{itr} = OptimalPeriod;
    Traditional{itr} = train;
    clear seg;
    
    %Joint Dist. Params. Estimation
    exp_pd_GCS = fitdist(OptimalPeriod(:,1),'Exp');
    exp_cdf_GCS = cdf(exp_pd_GCS,OptimalPeriod(:,1));
    gam_pd_GCS = fitdist(OptimalPeriod(:,2),'Gamma');
    gam_cdf_GCS = cdf(gam_pd_GCS,OptimalPeriod(:,2));
    paramhat_GCS = copulafit(family,[exp_cdf_GCS gam_cdf_GCS]);

    %Test LL
    loglikelihood_GCS(itr) = log(prod(copulapdf(family,[cdf(exp_pd_GCS,test(:,1)),cdf(gam_pd_GCS,test(:,2))],paramhat_GCS),'All'));
end
times(i) = cputime - start;
end

function loglikelihood = LL(x)
    global family lambda
    %marginal distribution fitting
    exp_pd = fitdist(x(:,1),'Exp');
    exp_cdf = cdf(exp_pd,x(:,1));
    exp_var = var(x(:,1));
    gam_pd = fitdist(x(:,2),'Gamma');
    gam_cdf = cdf(gam_pd,x(:,2));
    gam_var = var(x(:,2));

    %copula fitting
    paramhat = copulafit(family,[exp_cdf gam_cdf]);

    %loglikelihood
    loglikelihood = log(prod(copulapdf(family,[cdf(exp_pd,x(:,1)),cdf(gam_pd,x(:,2))],paramhat),'All')) - lambda / (exp_var + gam_var);
end

