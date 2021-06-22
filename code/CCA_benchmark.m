clc; clear all; close all;

%define copula family
global family lambda 
family = 'Clayton';
lambda = 100; %regularization parameter
num_rz = 10;
enforced_run = true; %continue the algorithm even though regularization is not good enough

result = zeros(num_rz,6);
result_GCS = zeros(num_rz,6);

for rz = 1:num_rz
    disp(rz);
    n = [300;300;400];
    alpha = [1;10;50];
    shape = [10;40;100]; 
    scale = [0.5;0.25;0.15];
    mu = [0.6628;1.0849;1.3785]; % [2,3,4]
    sigma = [0.2462;0.1655;0.1245]; % [0.5,0.5,0.5]
    
    data = benchmark_generator(n,alpha,mu,sigma,shape,scale);

    %define variables
    period = 1000;
    base = 400;
    cycle = 100;
    loglikelihood_GCS = zeros(6,1);
    loglikelihood = zeros(6,1);

        for itr=1:6
            %load data     
            train = data(1:base+cycle*(itr-1),:);
            test = data(base+cycle*(itr-1)+1:base+cycle*itr,:);

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
            OptimalPeriods{rz,itr} = OptimalPeriod;
            Traditional{rz,itr} = train;
            clear seg;

            %Joint Dist. Params. Estimation
            gam_pd_GCS = fitdist(OptimalPeriod(:,1),'Gamma');
            gam_cdf_GCS = cdf(gam_pd_GCS,OptimalPeriod(:,1));
            ln_pd_GCS = fitdist(OptimalPeriod(:,2),'Lognormal');
            ln_cdf_GCS = cdf(ln_pd_GCS,OptimalPeriod(:,2));
            paramhat_GCS = copulafit(family,[gam_cdf_GCS ln_cdf_GCS]);

            gam_pd = fitdist(train(:,1),'Gamma');
            gam_cdf = cdf(gam_pd,train(:,1));
            ln_pd = fitdist(train(:,2),'Lognormal');
            ln_cdf = cdf(ln_pd,train(:,2));
            paramhat = copulafit(family,[gam_cdf ln_cdf]);
            
            gam_pd_test = fitdist(test(:,1),'Gamma');
            gam_cdf_test = cdf(gam_pd_test,test(:,1));
            ln_pd_test = fitdist(test(:,2),'Lognormal');
            ln_cdf_test = cdf(ln_pd_test,test(:,2));
            paramhat_test = copulafit(family,[gam_cdf_test ln_cdf_test]);

            %Test LL
            loglikelihood_GCS(itr) = sum(log(copulapdf(family,[cdf(gam_pd_GCS,test(:,1)),cdf(ln_pd_GCS,test(:,2))],paramhat_GCS)));
            loglikelihood(itr) = sum(log(copulapdf(family,[cdf(gam_pd,test(:,1)),cdf(ln_pd,test(:,2))],paramhat)));
        end
    if (sum(loglikelihood_GCS < 0)  > 0) & (~enforced_run)
        break
    end
    result_GCS(rz,:) = loglikelihood_GCS;
    result(rz,:) = loglikelihood;
end

csvwrite(strcat('../outputs/bm_result_l_',num2str(lambda),'.csv'),result);
csvwrite(strcat('../outputs/bm_result_GCS_l_',num2str(lambda),'.csv'),result_GCS);

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

function x = benchmark_generator(n,alpha,mu,sigma,shape,scale)
    global family
    for i=1:length(n)
        u = copularnd(family,alpha(i),n(i));
        pd1 = makedist('Gamma','a',shape(i),'b',scale(i));
        pd2 = makedist('Lognormal','mu',mu(i),'sigma',sigma(i));
        if exist('x1','var')
            x1 = [x1;icdf(pd1,u(:,1))];
            x2 = [x2;icdf(pd2,u(:,2))];
        else
            x1 = icdf(pd1,u(:,1));
            x2 = icdf(pd2,u(:,2));
        end
    end
    x = [x1,x2];
end