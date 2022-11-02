clc; clear all; close all;

global copula_family var1_family var2_family
copula_family = 'Clayton';
var1_family = 'Gamma';
var2_family = 'Lognormal';

train = csvread(strcat('../../data/mhw/mhw_2017_train.csv'),1,0);
optimal = train(8:52,:);
test = csvread(strcat('../../data/mhw/mhw_2017_test.csv'),1,0);
traintest = cat(1,train,test);
optimaltest = cat(1,optimal,test);

plotter(train);
plotter(optimal);
plotter(test);
%plotter(traintest);
%plotter(optimaltest);

function plotter(x)
    global copula_family var1_family var2_family
    MaxVal = 1000;
    options = optimset('MaxFunEvals',MaxVal,'MaxIter',MaxVal,'Display','none');
    %marginal distribution fitting
    pd1 = @(u)gamlike([u(1),u(2)],(x(:,1)-u(3)));
    params1 = fminsearch(pd1,[1,1,0],options);
    cdf1 = cdf(makedist(var1_family,"a",params1(1),"b",params1(2)),x(:,1)-params1(3));   

    pd2 = @(u)lognlike([u(1),u(2)],(x(:,2)-u(3)));
    params2 = fminsearch(pd2,[1,1,0],options);
    cdf2 = cdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),x(:,2)-params2(3));   
    
    %copula fitting
    paramhat = copulafit(copula_family,[cdf1 cdf2]);
    
    %plot
    figure()
    u1 = linspace(0.001,0.999,1000);
    u2 = linspace(0.001,0.999,1000);
    x1 = icdf(makedist(var1_family,"a",params1(1),"b",params1(2)),u1)+params1(3); 
    x2 = icdf(makedist(var2_family,"mu",params2(1),"sigma",params2(2)),u2)+params2(3); 
    [U1,U2] = meshgrid(u1,u2);
    y = copulapdf(copula_family,[U1(:),U2(:)],0.0376);
    contour(x1,x2,reshape(y,size(U1)),[0.93,0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04]);
    hold on
    xlabel('duration');
    ylabel('intensity_{max}');
    grid on;
    cb = colorbar;
    cb.Label.String = 'Joint Density';
    set(gca,'xlim',[4,100],'ylim',[1.5,7],'xscale','log','yscale','log');
    legend('','test-ldown','test-lup','test-lbigup','');
    title('Traditional-CCA');
    %print(gcf,strcat('../../plots/mhw/mhwCopulaTrad.png'),'-dpng','-r300');
    %set(gca,'xlim',[4.99995,5.0002],'ylim',[2.2115662,2.21156645],'xscale','log','yscale','log');
    %print(gcf,strcat('../../plots/mhw/mhwCopulaTrad_zoom.png'),'-dpng','-r300');
end