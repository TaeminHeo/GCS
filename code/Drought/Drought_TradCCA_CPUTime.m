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
loglikelihood = zeros(fix((total-base)/cycle),1);

start = cputime;
for itr=1:fix((total-base)/cycle)
    %load data     
    train = csvread(strcat('../../data/b',num2str(base),'u',num2str(cycle),'_',num2str(itr-1),'.csv'),1,0);
    test = csvread(strcat('../../data/b',num2str(base),'u',num2str(cycle),'_',num2str(itr-1),'_test.csv'),1,0);

    %Result
    Traditional{itr} = train;
    
    %Joint Dist. Params. Estimation
    exp_pd = fitdist(train(:,1),'Exp');
    exp_cdf = cdf(exp_pd,train(:,1));
    gam_pd = fitdist(train(:,2),'Gamma');
    gam_cdf = cdf(gam_pd,train(:,2));
    paramhat = copulafit(family,[exp_cdf gam_cdf]);

    %Test LL
    loglikelihood(itr) = log(prod(copulapdf(family,[cdf(exp_pd,test(:,1)),cdf(gam_pd,test(:,2))],paramhat),'All'));
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

