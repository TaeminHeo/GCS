clc; close all;

f_marginal_pdf_compare = 1;
f_copula_compare = 1;

exp_pd_train = fitdist(train(:,1),'Exp')
mean(exp_pd_train)
var(exp_pd_train)
exp_cdf_train = cdf(exp_pd_train,train(:,1));
gam_pd_train = fitdist(train(:,2),'Gamma')
mean(gam_pd_train)
var(gam_pd_train)
gam_cdf_train = cdf(gam_pd_train,train(:,2));
paramhat_train = copulafit(family,[exp_cdf_train gam_cdf_train])

exp_pd_optimal = fitdist(OptimalPeriod(:,1),'Exp')
mean(exp_pd_optimal)
var(exp_pd_optimal)
exp_cdf_optimal = cdf(exp_pd_optimal,OptimalPeriod(:,1));
gam_pd_optimal = fitdist(OptimalPeriod(:,2),'Gamma')
mean(gam_pd_optimal)
var(gam_pd_optimal)
gam_cdf_optimal = cdf(gam_pd_optimal,OptimalPeriod(:,2));
paramhat_optimal = copulafit(family,[exp_cdf_optimal gam_cdf_optimal])

exp_pd_test = fitdist(test(:,1),'Exp')
mean(exp_pd_test)
var(exp_pd_test)
exp_cdf_test = cdf(exp_pd_test,test(:,1));
gam_pd_test = fitdist(test(:,2),'Gamma')
mean(gam_pd_test)
var(gam_pd_test)
gam_cdf_test = cdf(gam_pd_test,test(:,2));
paramhat_test = copulafit(family,[exp_cdf_test gam_cdf_test])


if f_marginal_pdf_compare == 1
    figure;
    hold on
    x = linspace(0,20,50);
    plot(x,cdf(exp_pd_train,x),'LineWidth',1);
    plot(x,cdf(exp_pd_optimal,x),'LineWidth',1);
    plot(x,cdf(exp_pd_test,x),'LineWidth',1);
    xlabel('Duration');
    ylabel('Probability');
    set(gca,'Fontsize',16);
    hold off
    
    axes('position',[.436 .24 .45 .55])
    set(gca,'Fontsize',12);
    box on
    hold on
    zs = 20;
    ze = 50;
    tmp = cdf(exp_pd_train,x);
    plot(x(zs:ze),tmp(zs:ze));
    tmp = cdf(exp_pd_optimal,x);
    plot(x(zs:ze),tmp(zs:ze));
    tmp = cdf(exp_pd_test,x);
    plot(x(zs:ze),tmp(zs:ze));
    axis tight
    legend('train set','optimal period','test set','Location','southeast');
    saveas(gcf,'../plots/GCS_marginal_pdf_comparison_duration.png');
      
    figure;
    hold on
    x = linspace(0,15,50);
    plot(x,cdf(gam_pd_train,x),'LineWidth',1);
    plot(x,cdf(gam_pd_optimal,x),'LineWidth',1);
    plot(x,cdf(gam_pd_test,x),'LineWidth',1);
    xlabel('Severity');
    ylabel('Probability');
    set(gca,'Fontsize',16);
    hold off
    
    axes('position',[.436 .24 .45 .55])
    set(gca,'Fontsize',12);
    box on
    hold on
    zs = 20;
    ze = 50;
    tmp = cdf(gam_pd_train,x);
    plot(x(zs:ze),tmp(zs:ze));
    tmp = cdf(gam_pd_optimal,x);
    plot(x(zs:ze),tmp(zs:ze));
    tmp = cdf(gam_pd_test,x);
    plot(x(zs:ze),tmp(zs:ze));
    axis tight
    legend('train set','optimal period','test set','Location','southeast');
    saveas(gcf,'../plots/GCS_marginal_pdf_comparison_severity.png');
end

if f_copula_compare == 1  
    lins_d = linspace(0,20,100);
    lins_s = linspace(0,15,100);
    [d,s] = meshgrid(lins_d,lins_s);
    [d_train,s_train] = meshgrid(cdf(exp_pd_train,lins_d),cdf(gam_pd_train,lins_s));
    [d_optimal,s_optimal] = meshgrid(cdf(exp_pd_optimal,lins_d),cdf(gam_pd_optimal,lins_s));
    [d_test,s_test] = meshgrid(cdf(exp_pd_test,lins_d),cdf(gam_pd_test,lins_s));
    yy_train = copulacdf(family,[d_train(:),s_train(:)],paramhat_train);
    yy_optimal = copulacdf(family,[d_optimal(:),s_optimal(:)],paramhat_optimal);
    yy_test = copulacdf(family,[d_test(:),s_test(:)],paramhat_test);
    
    figure;
    surf(d,s,reshape(yy_test-yy_train,100,100));
    xlabel('Duration');
    ylabel('Severity');
    title({'Probability Difference','(test set - train set)'});
    view(0,90);
    colorbar();
    caxis([-0.0114 0.0383]);
    saveas(gcf,'../plots/GCS_prob_compare_test_train.png');
    
    figure;
    surf(d,s,reshape(yy_test-yy_optimal,100,100));
    xlabel('Duration');
    ylabel('Severity');
    title({'Probability Difference','(test set - optimal period)'});
    view(0,90);
    colorbar();
    caxis([-0.0114 0.0383]);
    saveas(gcf,'../plots/GCS_prob_compare_test_opt.png');
    
    figure;
    surf(d,s,reshape(yy_test-yy_train,100,100));
    xlabel('Duration');
    ylabel('Severity');
    title({'Probability Difference','(test set - train set)'});
    view(-40,40);
    colorbar();
    caxis([-0.0114 0.0383]);
    zlim([-0.0114 0.0383]);
    saveas(gcf,'../plots/GCS_prob_compare2_test_train.png');
    
    figure;
    surf(d,s,reshape(yy_test-yy_optimal,100,100));
    xlabel('Duration');
    ylabel('Severity');
    title({'Probability Difference','(test set - optimal period)'});
    view(-40,40);
    colorbar();
    caxis([-0.0114 0.0383]);
    zlim([-0.0114 0.0383]);
    saveas(gcf,'../plots/GCS_prob_compare2_test_opt.png');
    
    figure;
    surf(d,s,reshape(yy_test-yy_train,100,100));
    hold on
    surf(d,d*0.5,d*0+1,'FaceColor','k','EdgeColor','k','LineWidth',2);
    surf(d,d*1.0,d*0+1,'FaceColor','k','EdgeColor','k','LineWidth',2);
    surf(d,d*1.5,d*0+1,'FaceColor','k','EdgeColor','k','LineWidth',2);
    caxis([-0.0114 0.0383]);
    view(0,90);
    xlim([0,20]);
    ylim([0,15]);
    xlabel('Duration');
    ylabel('Severity');
    colorbar();
    title({'Probability Difference','(test set - train set)'});
    saveas(gcf,'../plots/Combined_Average_Drought_Intensity_Domains_test_train.png');
    
    figure;
    surf(d,s,reshape(yy_test-yy_optimal,100,100));
    hold on
    surf(d,d*0.5,d*0+1,'FaceColor','k','EdgeColor','k','LineWidth',2);
    surf(d,d*1.0,d*0+1,'FaceColor','k','EdgeColor','k','LineWidth',2);
    surf(d,d*1.5,d*0+1,'FaceColor','k','EdgeColor','k','LineWidth',2);
    caxis([-0.0114 0.0383]);
    view(0,90);
    xlim([0,20]);
    ylim([0,15]);
    xlabel('Duration');
    ylabel('Severity');
    colorbar();
    title({'Probability Difference','(test set - optimal period)'});
    saveas(gcf,'../plots/Combined_Average_Drought_Intensity_Domains_test_opt.png');
   
    figure;
    scatter(test(:,1),test(:,2),'k');
    hold on
    xlim([0,20]);
    ylim([0,15]);
    xlabel('Duration');
    ylabel('Severity');
    title({'Average Drought Intensity Domains',''});
    dim = [0.8 0.3 0.1 0.1];
    str = 'light';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',12);
    dim = [0.75 0.8 0.1 0.1];
    str = 'moderate';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',12);
    dim = [0.53 0.8 0.1 0.1];
    str = 'severe';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',12);
    dim = [0.2 0.8 0.1 0.1];
    str = 'extreme';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',12);
    plot(lins_d,lins_d*0.5,'k');
    plot(lins_d,lins_d*1.0,'k');
    plot(lins_d,lins_d*1.5,'k');
    saveas(gcf,'../plots/Average_Drought_Intensity_Domains.png');
end




