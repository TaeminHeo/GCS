close all;

family = 'Clayton';
alpha = [1;10;50];
% gamma mean=[3.5,5,6.45] var=[2.45,2.5,2.7735]
shape = [10;40;100]; 
scale = [0.5;0.25;0.15];
% lognormal mean=[2,3,4] var=[0.5,0.5,0.5]
mu = [0.6628;1.0849;1.3785]; % [2,3,4]
sigma = [0.2462;0.1655;0.1245]; % [0.5,0.5,0.5]

figure()
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 800, 200]);
tiledlayout(1,3);
for i = 1:3
    nexttile
    pd1 = makedist('Gamma','a',shape(i),'b',scale(i));
    pd2 = makedist('Lognormal','mu',mu(i),'sigma',sigma(i));
    
    x1 = linspace(icdf(pd1,0.01),icdf(pd1,0.99),20);
    x2 = linspace(icdf(pd2,0.01),icdf(pd2,0.99),20);
    u1 = cdf(pd1,x1);
    u2 = cdf(pd2,x2);
    [U1,U2] = meshgrid(u1,u2);
    y = log(copulapdf(family,[U1(:),U2(:)],alpha(i)));
    contourf(x1,x2,reshape(y,20,20),20);
    xlabel('x_1');
    ylabel('x_2');
    grid on;
    cb = colorbar('southoutside');
end
%print(gcf,strcat('../../plots/Copulas.png'),'-dpng','-r300');