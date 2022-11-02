clc; clear all; close all;

%Get file names
itr_max = 4;
base = 20;
cycle = 10;

%define copula family
global family lambda 
family = 'Frank';
lambdas = [0.01];

for i=1:length(lambdas)
lambda = lambdas(i); %regularization parameter
do_plot = false;

%define variables
loglikelihood_GCS = zeros(itr_max,1);
loglikelihood = zeros(itr_max,1);

for itr=1:itr_max
    %load data     
    train = csvread(strcat('../../data/SalinitySST/b',num2str(base),'u',num2str(cycle),'_',num2str(itr-1),'.csv'),1,0);
    test = csvread(strcat('../../data/SalinitySST/b',num2str(base),'u',num2str(cycle),'_',num2str(itr-1),'_test.csv'),1,0);
    
    %greedy copula segmentation
    K = 1;
    seg{1} = train;
    while K > 0
        %K
        LL_all = 0;
    for j=1:K
        n = length(seg{j});
        LL_segorig = LL(seg{j});
        for k=4:n-4
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
    var1_pd_GCS = fitdist(OptimalPeriod(:,1),'Normal');
    var1_cdf_GCS = cdf(var1_pd_GCS,OptimalPeriod(:,1));
    var2_pd_GCS = fitdist(OptimalPeriod(:,2),'Normal');
    var2_cdf_GCS = cdf(var2_pd_GCS,OptimalPeriod(:,2));
    paramhat_GCS = copulafit(family,[var1_cdf_GCS var2_cdf_GCS]);

    var1_pd = fitdist(train(:,1),'Normal');
    var1_cdf = cdf(var1_pd,train(:,1));
    var2_pd = fitdist(train(:,2),'Normal');
    var2_cdf = cdf(var2_pd,train(:,2));
    paramhat = copulafit(family,[var1_cdf var2_cdf]);

    %Test LL
    loglikelihood_GCS(itr) = log(prod(copulapdf(family,[cdf(var1_pd_GCS,test(:,1)),cdf(var2_pd_GCS,test(:,2))],paramhat_GCS),'All'));
    loglikelihood(itr) = log(prod(copulapdf(family,[cdf(var1_pd,test(:,1)),cdf(var2_pd,test(:,2))],paramhat),'All'));
end

csvwrite(strcat('../../outputs/SalinitySST/ss_b',num2str(base),'u',num2str(cycle),'_result_l_',num2str(lambda),'.csv'),loglikelihood);
csvwrite(strcat('../../outputs/SalinitySST/ss_b',num2str(base),'u',num2str(cycle),'_result_GCS_l_',num2str(lambda),'.csv'),loglikelihood_GCS);

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
    var1_pd = fitdist(x(:,1),'Normal');
    var1_cdf = cdf(var1_pd,x(:,1));
    var1_var = var(x(:,1));
    var2_pd = fitdist(x(:,2),'Normal');
    var2_cdf = cdf(var2_pd,x(:,2));
    var2_var = var(x(:,2));

    %copula fitting
    paramhat = copulafit(family,[var1_cdf var2_cdf]);

    %loglikelihood
    loglikelihood = log(prod(copulapdf(family,[cdf(var1_pd,x(:,1)),cdf(var2_pd,x(:,2))],paramhat),'All')) - lambda / (var1_var + var2_var);
end

