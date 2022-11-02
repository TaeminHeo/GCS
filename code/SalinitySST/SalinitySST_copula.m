close all;

family = 'Gaussian';
%alpha = [10.1626;3.4887;6.8768]; % Frank
%alpha = [2.7705;1.3863;2.0156]; % Gumbel
alpha = [0.8819;0.5336;0.7588]; % Gaussian
% normal SC
mu1 = [-0.00670281;-0.0117963;-0.00316565]; 
sigma1 = [0.00513741;0.00177572;0.00346084];
% normal SSTA
mu2 = [-0.0280269;-0.20793;0.096906]; 
sigma2 = [0.201262;0.105773;0.150887];

figure()
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 800, 200]);
tiledlayout(1,3);
for i = 1:3
    nexttile
    pd1 = makedist('Normal','mu',mu1(i),'sigma',sigma1(i));
    pd2 = makedist('Normal','mu',mu2(i),'sigma',sigma2(i));
    
    %x1 = linspace(icdf(pd1,0.001),icdf(pd1,0.999),30);
    %x2 = linspace(icdf(pd2,0.001),icdf(pd2,0.999),30);
    %u1 = cdf(pd1,x1);
    %u2 = cdf(pd2,x2);
    u1 = linspace(0.01,0.99,30);
    u2 = linspace(0.01,0.99,30);
    x1 = icdf(pd1,u1);
    x2 = icdf(pd2,u2);
    [U1,U2] = meshgrid(u1,u2);
    y = copulapdf(family,[U1(:),U2(:)],alpha(i));
    contourf(x1,x2,log(reshape(y,30,30)),[-3,-2,-1,0,1,2,3],'ShowText','on');
    if i == 1
        hold on
        %plot(train(:,1),train(:,2),'w.');
    elseif i == 2
        hold on
        %plot(train(1:300,1),train(1:300,2),'w.');
        %plot(train(300:732,1),train(300:732,2),'k.');
    else
        hold on
        %plot(train(1:300,1),train(1:300,2),'k.');
        %plot(train(300:732,1),train(300:732,2),'w.');
    end
    xlabel('x_1');
    ylabel('x_2');
    grid on;
    xlim([-20*10^-3,5*10^-3]);
    ylim([-0.5,0.5]);
    cb = colorbar;
end
cb.Label.String = 'Log Joint Density';
