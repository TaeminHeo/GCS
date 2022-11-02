clc

train = csvread('../../data/Benchmark/toyproblem.csv');


% gamma mean=[5,10,15] sig=[2.5,2.5,2.5]
shape = [10;40;100]; 
scale = [0.5;0.25;0.15];

% lognormal mean=[2,3,4] sig=[0.5,0.5,0.5]
mu = [0.6628;1.0849;1.3785]; % [2,3,4]
sigma = [0.2462;0.1655;0.1245]; % [0.5,0.5,0.5]

st = [1,301,601];
ed = [300,600,1000];

for i=1:3
    pd1 = makedist('Gamma','a',shape(i),'b',scale(i));
    pd2 = makedist('Lognormal','mu',mu(i),'sigma',sigma(i));
    cdf1 = cdf(pd1,train(st(i):ed(i),1));
    cdf2 = cdf(pd2,train(st(i):ed(i),2));
    paramhat = copulafit('Clayton',[cdf1 cdf2])
end