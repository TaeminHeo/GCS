close all;

figure()
set(gca,'FontSize',12)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2)-200, 800, 450]);
subplot(2,2,1)
scatter(train(1:481,1),train(1:481,2),5,'black');
xlabel('x_1');
ylabel('x_2');
title('X_1 breakpoint at k=481');
grid on;
xlim([0,30])
ylim([0,6])

subplot(2,2,2)
scatter(train(482:722,1),train(482:722,2),5,'blue');
hold on
scatter(train(723:1000,1),train(723:1000,2),5,'black');
xlabel('x_1');
ylabel('x_2');
title('X_2 breakpoint at k=481');
grid on;
xlim([0,30])
ylim([0,6])

subplot(2,2,3)
scatter(train(1:481,1),train(1:481,2),5,'black');
hold on
scatter(train(482:722,1),train(482:722,2),5,'blue');
xlabel('x_1');
ylabel('x_2');
title('X_1 breakpoint at k=722');
grid on;
xlim([0,30])
ylim([0,6])

subplot(2,2,4)
scatter(train(723:1000,1),train(723:1000,2),5,'black');
xlabel('x_1');
ylabel('x_2');
title('X_2 breakpoint at k=722');
grid on;
xlim([0,30])
ylim([0,6])
%print(gcf,strcat('../plots/greedy_compare.png'),'-dpng','-r300');

