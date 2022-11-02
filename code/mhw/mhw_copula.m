close all;

figure()
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
tiledlayout(1,2);

x1 = linspace(5,100,30);
x2 = linspace(1,6,30);

levels = [0.2,0.4,0.6,0.8,1,1.2,1.4];

nexttile
u1 = cdf(var1_pd,x1);
u2 = cdf(var2_pd,x2);
[U1,U2] = meshgrid(u1,u2);
y = copulapdf(copula_family,[U1(:),U2(:)],paramhat);
contourf(x1,x2,reshape(y,size(U1)),levels);
hold on
plot(test(:,1),test(:,2),'k.');
plot(train(:,1),train(:,2),'r.');
xlabel('var_1');
ylabel('var_2');
grid on;
cb = colorbar;

nexttile
u1 = cdf(var1_pd_GCS,x1);
u2 = cdf(var2_pd_GCS,x2);
[U1,U2] = meshgrid(u1,u2);
y = copulapdf(copula_family,[U1(:),U2(:)],paramhat_GCS);
contourf(x1,x2,reshape(y,size(U1)),levels);
hold on
plot(test(:,1),test(:,2),'k.');
plot(OptimalPeriod(:,1),OptimalPeriod(:,2),'r.');
xlabel('var_1');
ylabel('var_2');
grid on;
cb = colorbar;

cb.Label.String = 'Joint Density';
