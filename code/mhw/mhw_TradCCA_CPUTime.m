clc; clear all; close all;

%define copula family
global copula_family var1_family var2_family lambda 
copula_family = 'Clayton';
var1_family = 'Gamma';
var2_family = 'Lognormal';
times = zeros(1,30);

%load data     
train = csvread(strcat('../../data/mhw/mhw_2017_train.csv'),1,0);
test = csvread(strcat('../../data/mhw/mhw_2017_test.csv'),1,0);

for i=1:30
    start = cputime;
    
    %Joint Dist. Params. Estimation
    pd1 = @(u)gamlike([u(1),u(2)],(train(:,1)-u(3)));
    params1 = fminsearch(pd1,[1,1,0]);
    cdf1 = cdf(makedist(var1_family,"a",params1(1),"b",params1(2)),train(:,1)-params1(3));   
    pd2 = @(u)lognlike([u(1),u(2)],(train(:,2)-u(3)));
    params2 = fminsearch(pd2,[1,1,0]);
    cdf2 = cdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),train(:,2)-params2(3));   
    paramhat = copulafit(copula_family,[cdf1 cdf2]);

    times(i) = cputime - start;
end

mean(times)
std(times)
