close all;
train = csvread('../data/toyproblem.csv');

x = train;
family = 'Clayton';
alpha = [1;10;50];
rate = [5;10;15];
mu = [0.6628;1.0849;1.3785]; % [2,3,4]
sigma = [0.2462;0.1655;0.1245]; % [0.5,0.5,0.5]

figure()
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 800, 200]);
tiledlayout(1,3);
for i = 1:3
    nexttile
    pd1 = makedist('Poisson','lambda',rate(i));
    pd2 = makedist('Lognormal','mu',mu(i),'sigma',sigma(i));
    
    x1 = linspace(icdf(pd1,0.1),icdf(pd1,0.9),10);
    x2 = linspace(icdf(pd2,0.1),icdf(pd2,0.9),10);
    u1 = cdf(pd1,x1);
    u2 = cdf(pd2,x2);
    [U1,U2] = meshgrid(u1,u2);
    y = copulapdf('Clayton',[U1(:),U2(:)],alpha(i));
    contourf(x1,x2,reshape(y,10,10));
    xlabel('x_1');
    ylabel('x_2');
    grid on;
    cb = colorbar;
end
cb.Label.String = 'Joint Density';
print(gcf,strcat('../plots/Copulas.png'),'-dpng','-r300');