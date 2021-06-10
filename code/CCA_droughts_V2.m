clc; clear all; close all;

%Get file names
total = 120;
base = 20;
cycle = 10;
fnames = dir(strcat('../data/b',num2str(base),'u',num2str(cycle),'_*_V2.csv'));
fnames = {fnames.name};

%define copula family
global family lambda 
family = 'Gumbel';
lambda = 100; %150
do_plot = true;

%define variables
loglikelihood_GCS = zeros(fix((total-base)/cycle),1);
loglikelihood = zeros(fix((total-base)/cycle),1);

for itr=1:fix((total-base)/cycle)
    %load data
    train = csvread(strcat('../data/',string(fnames(2*(itr-1)+1))),1,0);
    test = csvread(strcat('../data/',string(fnames(2*itr))),1,0);

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
    ln_pd_GCS = fitdist(OptimalPeriod(:,3),'Gamma');
    ln_cdf_GCS = cdf(ln_pd_GCS,OptimalPeriod(:,3));
    paramhat_GCS = copulafit(family,[exp_cdf_GCS ln_cdf_GCS]);

    exp_pd = fitdist(train(:,1),'Exp');
    exp_cdf = cdf(exp_pd,train(:,1));
    gam_pd = fitdist(train(:,2),'Gamma');
    gam_cdf = cdf(gam_pd,train(:,2));
    ln_pd = fitdist(train(:,3),'Gamma');
    ln_cdf = cdf(ln_pd,train(:,3));
    paramhat = copulafit(family,[exp_cdf ln_cdf]);

    %Test LL
    loglikelihood_GCS(itr) = log(prod(copulapdf(family,[cdf(exp_pd_GCS,test(:,1)),cdf(ln_pd_GCS,test(:,3))],paramhat_GCS),'All'));
    loglikelihood(itr) = log(prod(copulapdf(family,[cdf(exp_pd,test(:,1)),cdf(ln_pd,test(:,3))],paramhat),'All'));
end

csvwrite(strcat('../outputs/dr_b',num2str(base),'u',num2str(cycle),'_result_l_',num2str(lambda),'_V2.csv'),loglikelihood);
csvwrite(strcat('../outputs/dr_b',num2str(base),'u',num2str(cycle),'_result_GCS_l_',num2str(lambda),'_V2.csv'),loglikelihood_GCS);

if do_plot
    bar((0:(fix((total-base)/cycle)-1))*10,100*(loglikelihood_GCS - loglikelihood) ./ abs(loglikelihood_GCS),'Facecolor',[0.9,0.9,0.9]);
    grid on
    xlim([-10,100]);
    xlabel('Number of additional years following the first 20 years');
    ylabel('\delta_L_L (%)');
    print(gcf,strcat('../plots/DroughtCCAResults_',num2str(lambda),'_V2.png'),'-dpng','-r300');
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
    ln_pd = fitdist(x(:,3),'Gamma');
    ln_cdf = cdf(ln_pd,x(:,3));
    ln_var = var(x(:,3));
    %copula fitting
    paramhat = copulafit(family,[exp_cdf ln_cdf]);

    %loglikelihood
    loglikelihood = log(prod(copulapdf(family,[cdf(exp_pd,x(:,1)),cdf(ln_pd,x(:,3))],paramhat),'All')) - lambda / (exp_var + ln_var);
end

