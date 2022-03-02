clc; clear all; close all;

global family
family = 'Clayton';

n = [300;300;400];
alpha = [1;10;50];
% gamma mean=[5,10,15] sig=[2.5,2.5,2.5]
shape = [10;40;100]; 
scale = [0.5;0.25;0.15];

% lognormal mean=[2,3,4] sig=[0.5,0.5,0.5]
mu = [0.6628;1.0849;1.3785]; % [2,3,4]
sigma = [0.2462;0.1655;0.1245]; % [0.5,0.5,0.5]

data = benchmark_generator(n,alpha,mu,sigma,shape,scale);
csvwrite('../../data/toyproblem.csv',data);

figure()
plot(data(:,1))

figure()
plot(data(:,2))

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