clc; clear all; close all;

%Get file names
total = 120;
base = 20;
cycle = 10;
fnames = dir(strcat('../../data/b',num2str(base),'u',num2str(cycle),'_*.csv'));
fnames = {fnames.name};

%define copula family
global family lambda 
family = 'Frank';
lambdas = [50,100,110,120,130,140,150,160,170,180,190,200,250];

for i=1:length(lambdas)
lambda = lambdas(i); %regularization parameter
do_plot = false;

%define variables
loglikelihood_GCS = zeros(fix((total-base)/cycle),1);
loglikelihood = zeros(fix((total-base)/cycle),1);

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
    norm_pd_GCS = fitdist(OptimalPeriod(:,1),'Normal');
    norm_cdf_GCS = cdf(norm_pd_GCS,OptimalPeriod(:,1));
    gam_pd_GCS = fitdist(OptimalPeriod(:,2),'Gamma');
    gam_cdf_GCS = cdf(gam_pd_GCS,OptimalPeriod(:,2));
    paramhat_GCS = copulafit(family,[norm_cdf_GCS gam_cdf_GCS]);

    norm_pd = fitdist(train(:,1),'Normal');
    norm_cdf = cdf(norm_pd,train(:,1));
    gam_pd = fitdist(train(:,2),'Gamma');
    gam_cdf = cdf(gam_pd,train(:,2));
    paramhat = copulafit(family,[norm_cdf gam_cdf]);

    %Test LL
    loglikelihood_GCS(itr) = log(prod(copulapdf(family,[cdf(norm_pd_GCS,test(:,1)),cdf(gam_pd_GCS,test(:,2))],paramhat_GCS),'All'));
    loglikelihood(itr) = log(prod(copulapdf(family,[cdf(norm_pd,test(:,1)),cdf(gam_pd,test(:,2))],paramhat),'All'));
end

csvwrite(strcat('../../outputs/dr_b',num2str(base),'u',num2str(cycle),'_result_l_',num2str(lambda),'.csv'),loglikelihood);
csvwrite(strcat('../../outputs/dr_b',num2str(base),'u',num2str(cycle),'_result_GCS_l_',num2str(lambda),'.csv'),loglikelihood_GCS);

if do_plot
    bar((0:(fix((total-base)/cycle)-1))*10,100*(loglikelihood_GCS - loglikelihood) ./ abs(loglikelihood_GCS),'Facecolor',[0.9,0.9,0.9]);
    grid on
    xlim([-10,100]);
    xlabel('Number of additional years following the first 20 years');
    ylabel('\delta_L_L (%)');
    %print(gcf,strcat('../../plots/DroughtCCAResults_',num2str(lambda),'.png'),'-dpng','-r300');
end
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

