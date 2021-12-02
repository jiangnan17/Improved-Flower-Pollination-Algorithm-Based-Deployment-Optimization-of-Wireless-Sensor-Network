

%三角形和菱形障碍物
geshu_x = [20,25,30,35,40,45,50];
init_y = [54.21,61.09,67.44,73.07,77.24,79.15,82.61]/100;
ga_y = [80.71,87.40,94.69,96.74,98.23,99.17,100]/100;
pso_y = [83.37,92.37,95.13,96.79,98.56,99.31,100]/100;
dea_y = [85.08,93.03,95.79,97.62,99.28,99.50,100]/100;
fa_y = [85.69,93.57,96.96,98.61,99.28,100,100]/100;
ifa_y = [86.02,94.25,97.12,98.73,99.34,100,100]/100;



figure(1);
% plot(geshu_x,init_y,'color','k');
% hold on;
plot(geshu_x,ga_y,'color','r');
hold on;
plot(geshu_x,pso_y,'color','g');
hold on;
plot(geshu_x,dea_y,'color','b');
hold on;
plot(geshu_x,fa_y,'color','c');
hold on;
plot(geshu_x,ifa_y,'color','m');
hold on;
legend('ga','pso','dea','fa','ifa');
hold on;



