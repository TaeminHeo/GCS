clc; clear all; close all;

%define copula family
global family lambda 
family = 'Clayton';
lambda = 500; %regularization parameter
num_rz = 10;
enforced_run = true; %continue the algorithm even though regularization is not good enough
do_plot = false;

result = zeros(num_rz,6);
result_GCS = zeros(num_rz,6);

for rz = 1:num_rz
    disp(rz);
    n = [400;300;300];
    alpha = [1;10;50];
    rate = [5;10;15];
    mu = [0.6628;1.0849;1.3785]; % [2,3,4]
    sigma = [0.2462;0.1655;0.1245]; % [0.5,0.5,0.5]
    data = benchmark_generator(n,alpha,rate,mu,sigma);

    %define variables
    period = 1000;
    base = 400;
    cycle = 100;
    loglikelihood_GCS = zeros(6,1);
    loglikelihood = zeros(6,1);

        for itr=1:6
            %load data     
            train = data(1:400+cycle*(itr-1),:);
            test = data(400+cycle*(itr-1)+1:400+cycle*itr,:);

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
            poi_pd_GCS = fitdist(OptimalPeriod(:,1),'Poisson');
            poi_cdf_GCS = cdf(poi_pd_GCS,OptimalPeriod(:,1));
            ln_pd_GCS = fitdist(OptimalPeriod(:,2),'Lognormal');
            ln_cdf_GCS = cdf(ln_pd_GCS,OptimalPeriod(:,2));
            paramhat_GCS = copulafit(family,[poi_cdf_GCS ln_cdf_GCS]);

            poi_pd = fitdist(train(:,1),'Poisson');
            poi_cdf = cdf(poi_pd,train(:,1));
            ln_pd = fitdist(train(:,2),'Lognormal');
            ln_cdf = cdf(ln_pd,train(:,2));
            paramhat = copulafit(family,[poi_cdf ln_cdf]);
            
            poi_pd_test = fitdist(test(:,1),'Poisson');
            poi_cdf_test = cdf(poi_pd_test,test(:,1));
            ln_pd_test = fitdist(test(:,2),'Lognormal');
            ln_cdf_test = cdf(ln_pd_test,test(:,2));
            paramhat_test = copulafit(family,[poi_cdf_test ln_cdf_test]);

            %Test LL
            loglikelihood_GCS(itr) = log(prod(copulapdf(family,[cdf(poi_pd_GCS,test(:,1)),cdf(ln_pd_GCS,test(:,2))],paramhat_GCS),'All'));
            loglikelihood(itr) = log(prod(copulapdf(family,[cdf(poi_pd,test(:,1)),cdf(ln_pd,test(:,2))],paramhat),'All'));
        end
    if (sum(loglikelihood_GCS < 0)  > 0) & (~enforced_run)
        break
    end
    result_GCS(rz,:) = loglikelihood_GCS;
    result(rz,:) = loglikelihood;
end

csvwrite(strcat('../outputs/bm_result_l_',num2str(lambda),'.csv'),result);
csvwrite(strcat('../outputs/bm_result_GCS_l_',num2str(lambda),'.csv'),result_GCS);

if do_plot
    diff = 100 * (result_GCS - result) ./ result_GCS;
    mean_diff = mean(diff);
    std_diff = std(diff);

    x=1:1:6;
    patch([x fliplr(x)], [min(diff) fliplr(max(diff))],[0.8 0.8 0.8]);
    hold on
    line1 = plot(mean_diff,'k');
    line2 = plot(min(diff),'k');
    line3 = plot(max(diff),'k');
    xticks([1 2 3 4 5 6]);
    xlabel('m');
    ylabel('\delta_L_L (%)');
    %saveas(gcf,'../plots/BenchmarkCCAResults.png');
end

function loglikelihood = LL(x)
    global family lambda
    %marginal distribution fitting
    poi_pd = fitdist(x(:,1),'Poisson');
    poi_cdf = cdf(poi_pd,x(:,1));
    poi_var = var(x(:,1));
    ln_pd = fitdist(x(:,2),'Lognormal');
    ln_cdf = cdf(ln_pd,x(:,2));
    ln_var = var(x(:,2));

    %copula fitting
    paramhat = copulafit(family,[poi_cdf ln_cdf]);

    %loglikelihood
    loglikelihood = log(prod(copulapdf(family,[cdf(poi_pd,x(:,1)),cdf(ln_pd,x(:,2))],paramhat),'All')) - lambda / (poi_var + ln_var);
end

function x = benchmark_generator(n,alpha,rate,mu,sigma)
    global family
    for i=1:length(n)
        u = copularnd(family,alpha(i),n(i));
        pd1 = makedist('Poisson','lambda',rate(i));
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