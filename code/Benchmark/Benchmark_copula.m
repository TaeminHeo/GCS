close all;

pd1 = fitdist(OptimalPeriod(:,1),'Gamma');
cdf1 = cdf(pd1,OptimalPeriod(:,1));
var1 = var(OptimalPeriod(:,1));
pd2 = fitdist(OptimalPeriod(:,2),'Lognormal');
cdf2 = cdf(pd2,OptimalPeriod(:,2));
var2 = var(OptimalPeriod(:,2));

%copula fitting
paramhat = copulafit(family,[cdf1 cdf2]);

figure()
set(gca,'FontSize',12)
pos = get(gcf, 'Position');   
x1 = linspace(icdf(pd1,0.1),icdf(pd1,0.9),20);
x2 = linspace(icdf(pd2,0.1),icdf(pd2,0.9),20);
u1 = cdf(pd1,x1);
u2 = cdf(pd2,x2);
[U1,U2] = meshgrid(u1,u2);
y = copulapdf(family,[U1(:),U2(:)],alpha(i));
contourf(x1,x2,reshape(y,20,20),20);
xlabel('x_1');
ylabel('x_2');
grid on;
cb = colorbar;
cb.Label.String = 'Joint Density';

