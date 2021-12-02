%%主程序
clc;
clear ;
close all;
%删除相应的文件

global N;
global M;
global L;
global W;
global Grid_cen_x;
global Grid_cen_y;
global Grid_cen_x_and_y;
global ger;

p=0.8;%判断是否是全局优化还是局部优化

L = 50;%长
W = 50;%宽
%假设1平方米一个网格
M = 2500;%网格总数
r_max = 7;%感知半径为5
r_mid = 6;
r_min = 5;
energy_max = 100;%最大的能量
energy_mid = 90;
energy_min = 80;

per_sersons_radius_type = [r_max,r_mid,r_min];
%假设大、中为5，剩下为小
N = 25;%30个传感器节点
sizepop = 50;%种群规模
dimension = 2;% 空间维数  前行放x、y，第三行放半径
ger = 10;% 最大迭代次数
pos_limit = [0, 50];            % 设置位置参数限制
%个数限制
r_max_num = 1;%序号为1-5
r_mid_num = 2;%序号为6-10
r_min_num = N - r_max_num - r_mid_num; %序号为11-N


struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);%结构体类型

struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);%临时的一个种群

energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);
%求出梯形的四个点
syms x y;%先定义一个变量
%左上角
k1 = 1;
b1 = 35;
x1_up = solve(k1*x+b1==50,x);%左上角的斜线的上个交点
y1_down = solve(k1*0+b1==y,y);

%左下角
k2 = -1;
b2 = 15;
y2_up = solve(k2*0+b2==y,y);
x2_down = solve(k2*x+b2==0,x);


%右上角
k3 = -1;
b3 = 85;
x3_up = solve(k3*x+b3==50,x);
y3_down = solve(k3*50+b3==y,y);

%右下角
k4 = 1;
b4 = -35;
y4_up = solve(k4*50+b4==y,y);
x4_down = solve(k4*x+b4==0,x);

%以下数据验证完毕，完全正确
point = zeros(8,2);%存储这些点  从左  从上往下
point(1,:) = [x1_up,50];
point(2,:) = [0,y1_down];
point(3,:) = [0,y2_up];
point(4,:) = [x2_down,0];
point(5,:) = [x3_up,50];
point(6,:) = [50,y3_down];
point(7,:) = [50,y4_up];
point(8,:) = [x4_down,0];

%菱形的计算
point_diamond = zeros(2,4);%菱形的四个点，方位是顺时针 第一列为上 二列为右
%求出新菱形形的四个点
syms x y;%先定义一个变量
%左上角
k5 = 1;
b5 = 10;
%别搞什么计算了  直接可以看出来

point_diamond(1,1) = 25;
point_diamond(2,1) = 35;

%右上角
k6 = -1;
b6 = 60;

point_diamond(1,2) = 35;
point_diamond(2,2) = 25;

%右下角
k7 = 1;
b7 = -10;

point_diamond(1,3) = 25;
point_diamond(2,3) = 15;

%左下角
k8 = 1;
b8 = 40;

point_diamond(1,4) = 15;
point_diamond(2,4) = 25;

load struct_pop_public.mat;%加载该种群
struct_pops = struct_pop_public;%得到种群数据

load struct_first_init_public.mat%加载最开始的一个个体数据
struct_first_init = struct_first_init_public;%得到初始化个体数据


%%初始的部署后画图  拿第一个粒子拿去初始画图

%求网格中心坐标
X_mat = (0:1:50);%x矩阵
Y_mat = (0:1:50);%y矩阵
Grid_cen_x = zeros(1,L/1);%网格中心点x坐标
Grid_cen_y = zeros(1,W/1);%网格中心点y坐标
%前后两者相加之和除以2
for i=1:L/1
    Grid_cen_x(i) = (X_mat(i)+X_mat(i+1))/2;
    Grid_cen_y(i) = (Y_mat(i)+Y_mat(i+1))/2;
end


%%把横纵坐标丢到一个二维矩阵当中
%用于转坐标  第一行放x轴 第二行放y轴，同一列放一个点坐标
%且先存靠近x坐标的第一行，然后往上存第二行
%网格中心坐标
Grid_cen_x_and_y = zeros(L,W,2);%共2500个网格中心，但是每个网格中心有x,y
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);%1代表x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);%把y坐标放到第二行
    end
end


x_pos = struct_first_init.per(1,:);%第一个粒子   粒子即是解  的x坐标
y_pos = struct_first_init.per(2,:);%第一个粒子   粒子即是解  的y坐标
sersors_r = struct_first_init.radius;%第一个粒子 




%进行画图
figure(1);
draw_circle(x_pos,y_pos,sersors_r);
title('初始化部署图');
hold on;





%把随机部署的点的坐标放到矩阵当中
sensor_mat = zeros(2,N);%预分配内存
for i=1:N
    sensor_mat(1,i) = x_pos(i);
    sensor_mat(2,i) = y_pos(i);
end

%存放三种节点的数目
per_sersons_num = struct_pops(1).sersons_num;
%存放三种节点的半径


%%初始化得到联合概率和节点浪费率
[cover_rate,waste_rate] =  get_Grid_cover_unit_and_rate_waste(sensor_mat,sersors_r,per_sersons_num,per_sersons_radius_type);
disp(['初始化的覆盖率：',num2str(cover_rate)]);
disp(['初始化的浪费率：',num2str(waste_rate)]);



%%计算连通性  第一个节点初始化时的连通性
is_connec = get_connection(sensor_mat,sersors_r);
if is_connec==1
    disp('连通');
else
    disp('非连通');
end




%%初始化种群历史值为  无穷小   inf为无穷大
best_fitness = -inf;                         % 种群历史最佳适应度  
struct_best_indivi = struct_pop_per;                 % 保存优秀个体
struct_best_indivi_fitness = struct('cover_rate',[],'waste_rate',[],'energy_rate',[],'function_rate',[]);%适应值结构体类型
struct_best_indivi_fitness_all = repmat(struct_best_indivi_fitness,[1 1]);%预分配内存
%%在上面已经画图了
%%以上为画出初始化时，50个粒子开始的位置

ifa_cover_fitness = zeros(sizepop,1);%个体覆盖率
ifa_waste_fitness = zeros(sizepop,1);%个体浪费率
ifa_energy_fitness = zeros(sizepop,1);%个体能耗率  化成1/xx形式 
ifa_function_fitness = zeros(sizepop,1);%三个函数的结合
weight_rate = [1,0,0];%设置权重，覆盖率应该是最大的


%% 群体更新
iter = 1;
record_ger = zeros(ger, 1);          % 记录每次迭代中最好的适应值 
record_pop_ave = zeros(ger,1);       % 记录种群适应值的平均值





while iter <= ger

    %计算该种群中每个个体的适应值

    %求各种适应值   记得我们这里求最大的  所以得化成求最大的
    for k=1:sizepop
        sensor_mat(1,:) = struct_pops(k).per(1,:);
        sensor_mat(2,:) = struct_pops(k).per(2,:);
        [ifa_cover_fitness(k,1), ifa_waste_fitness(k,1)] = get_Grid_cover_unit_and_rate_waste(sensor_mat,struct_pops(k).radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
        [~,ifa_energy_fitness(k,1)] = get_energy_consume(struct_first_init.per,struct_pops(k).per,struct_pops(k).radius,struct_pops(k).energy_init);
        ifa_function_fitness(k,1) = weight_rate(1,1) * ifa_cover_fitness(k,1) + weight_rate(1,2)*(1-ifa_energy_fitness(k,1)) + weight_rate(1,3) * (1 - ifa_waste_fitness(k,1));
    end
    
   
    
    [ifa_function_fitness_sort,order_index] = sort(ifa_function_fitness);%对适应值进行排序
    disp(ifa_function_fitness_sort);%用于临时的打印已经排好序号的适应值 便于查看
    

    %更新最优适应值和最优个体
    if ifa_function_fitness(order_index(sizepop,1),1) > best_fitness
        
        %用结构体保存相应的各种适应值  已经排好序了 所以是最后一个
        struct_best_indivi_fitness_all.cover_rate = ifa_cover_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.waste_rate = ifa_waste_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.energy_rate = ifa_energy_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.function_rate = ifa_function_fitness(order_index(sizepop,1),1);
        
        %保存最优秀的个体
        struct_best_indivi = struct_pops(order_index(sizepop,1));
        %保存最优适应值
        best_fitness = ifa_function_fitness(order_index(sizepop,1),1);
    end
    
    record_ger(iter,1) = best_fitness;%保存该代中最优的适应度值
    
    sum_fitness = sum(ifa_function_fitness);
    record_pop_ave(iter,1) = sum_fitness/sizepop;%算出该种群的平均适应度
    
    %用一个临时种群去操作
    struct_pops_new = struct_pops;

    
%   %进行特殊的处理
    num_swap = 1;
    for j=1:sizepop
        rand_index_swap_best = randperm(N,num_swap);
        %从最好的个体那里获取一部分节点
        for z=1:num_swap
            struct_pops_new(j).per(:,rand_index_swap_best(1,z)) = struct_best_indivi.per(:,rand_index_swap_best(1,z));
        end
    end
    
    
    
    b = 1-(1-(((ger-iter)/ger).^2)).^0.5;
    for i=1:sizepop 
        fitness_current = ifa_function_fitness(i,1);%得到当前个体的适应值
         if rand(1,1) < p               %局部搜索
            index_rand = randperm(sizepop,2);%随机选择两个
            epsilon = b*rand(2,N);%相当于是变化的步长吧
            %%得到新的中间体
            struct_temp_indivi_lo = struct_pops_new(i);
            index_N_serson = randperm(N,N);%所有结点进行处理
            for j=1:N
                struct_temp_indivi_lo.per(:,index_N_serson(1,j)) = struct_pops_new(i).per(:,index_N_serson(1,j)) + epsilon(:,j) .*(struct_pops_new(index_rand(1,1)).per(:,index_N_serson(1,j))- struct_pops_new(index_rand(1,2)).per(:,index_N_serson(1,j)));
                
                
                %同时进行越界处理  2表示x,y坐标
                for k=1:2
                     if struct_temp_indivi_lo.per(k,index_N_serson(1,j)) < pos_limit(1,1) || struct_temp_indivi_lo.per(k,index_N_serson(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_lo.per(k,index_N_serson(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %处理障碍物
                %左上
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %左下
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2表示相加
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2);
                     end
                 end


                 %处理第三个区域
                 %右上角
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(5,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(6,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3);
                     end
                 end

                 %处理第四区域
                 %右下角
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(8,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4);
                     end
                 end

                %处理菱形
                %画四个区域划分
                %处理左上斜边
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,4)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,4)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,1))
                     if (k5 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b5) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%相加%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k5 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %右上
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,1))
                     if (k6 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b6) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2表示相加
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k6 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b6);
                     end
                 end


                 %处理第三个区域
                 %右下角
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,3)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,3)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,2))
                     if (k7 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b7) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k7 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b7);
                     end
                 end

                 %处理第四区域
                 %左下角
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,4)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,3)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,3)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,4))
                     if (k8 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b8) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k8 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b8);
                     end
                 end
            end
            
            
            %用最好的个体进行下处理
            lo_swap_num = 1;
            lo_serson_index = randperm(N,lo_swap_num);
            for k=1:lo_swap_num
                struct_temp_indivi_lo.per(:,lo_serson_index(1,k)) = struct_best_indivi.per(:,lo_serson_index(1,k));
            end
           
            %局部搜索得到的适应值

             %得到新的适应值
            [temp_lo_cover_fitness, temp_lo_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_temp_indivi_lo.per,struct_temp_indivi_lo.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_lo_energy_fitness] = get_energy_consume(struct_first_init.per,struct_temp_indivi_lo.per,struct_temp_indivi_lo.radius,struct_temp_indivi_lo.energy_init);
            fitness_temp_function_lo = weight_rate(1,1) * temp_lo_cover_fitness + weight_rate(1,2)*(1-temp_lo_energy_fitness) + weight_rate(1,3) * (1 - temp_lo_waste_fitness);

            %%判断它的适应值是否超过了之前的  如果是则进行替换
            if fitness_temp_function_lo > fitness_current
                struct_pops_new(i) = struct_temp_indivi_lo;%进行替换
                
            end
         else%进行全局的搜索
            L1 = Levy2(dimension,iter,ger);%得到步长
            index_rand = randperm(N,N);%随机选择两个
            %%得到新的中间体
            struct_temp_indivi_glo = struct_pops_new(i);
            for j=1:N
                struct_temp_indivi_glo.per(:,index_rand(1,j)) = struct_pops_new(i).per(:,index_rand(1,j)) + L1' .*(struct_pops_new(i).per(:,index_rand(1,j))- struct_best_indivi.per(:,index_rand(1,j)));%得到临时花朵
            
                 %同时进行越界处理  2表示x,y坐标
                for k=1:2
                     if struct_temp_indivi_glo.per(k,index_rand(1,j)) < pos_limit(1,1) || struct_temp_indivi_glo.per(k,index_rand(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_glo.per(k,index_rand(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %处理障碍物
                %左上
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1) <= struct_temp_indivi_glo.per(2,index_rand(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %左下
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2表示相加
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2);
                     end
                 end


                 %处理第三个区域
                 %右上角
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(5,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(6,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3);
                     end
                 end

                 %处理第四区域
                 %右下角
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(8,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4);
                     end
                 end

                %处理菱形
                %画四个区域划分
                %处理左上斜边
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,4)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,4)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,1))
                     if (k5 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b5) >= struct_temp_indivi_glo.per(2,index_rand(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%相加%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k5 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %右上
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,1))
                     if (k6 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b6) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2表示相加
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k6 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b6);
                     end
                 end


                 %处理第三个区域
                 %右下角
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,3)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,3)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,2))
                     if (k7 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b7) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k7 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b7);
                     end
                 end

                 %处理第四区域
                 %左下角
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,4)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,3)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,3)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,4))
                     if (k8 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b8) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k8 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b8);
                     end
                 end
            
            end
            
            %用最好的个体进行下处理
            glo_swap_num = 2;
            glo_serson_index = randperm(N,glo_swap_num);
            for k=1:glo_swap_num
                struct_temp_indivi_glo.per(:,glo_serson_index(1,k)) = struct_best_indivi.per(:,glo_serson_index(1,k));
            end
            
             %得到新的适应值
            [temp_glo_cover_fitness, temp_glo_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_temp_indivi_glo.per,struct_temp_indivi_glo.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_glo_energy_fitness] = get_energy_consume(struct_first_init.per,struct_temp_indivi_glo.per,struct_temp_indivi_glo.radius,struct_temp_indivi_glo.energy_init);
            fitness_temp_function_glo = weight_rate(1,1) * temp_glo_cover_fitness + weight_rate(1,2)*(1-temp_glo_energy_fitness) + weight_rate(1,3) * (1 - temp_glo_waste_fitness);

            %%判断它的适应值是否超过了之前的  如果是则进行替换
            if fitness_temp_function_glo > fitness_current
                struct_pops_new(i) = struct_temp_indivi_glo;%进行替换
                
            end
           
         end
        
        %局部的交叉
        if rand(1,1) < 0.4;
            rand_swap_p = b*rand(1,1);%用于交换信息
            flower_rand_one_va = randperm(sizepop,2);
            %得到一个附近的flower进行交叉
            for kk=1:2
                if i ~= flower_rand_one_va(1,kk)
                    swap_index = flower_rand_one_va(1,kk);
                    break;
                end
            end
            %进行交叉
            %交换
            p_swap = b*rand(2,N);
            struct_flower_current_ac = struct_pops_new(i);
            struct_flower_one_ac = struct_pops_new(swap_index);
            ac_swap_num = 1;
            ac_swap_index = randperm(N,ac_swap_num);
            %交换
            for z=1:ac_swap_num
               
                struct_flower_current_ac.per(:,ac_swap_index(1,z)) = p_swap(:,ac_swap_index(1,z)) .* struct_pops_new(i).per(:,ac_swap_index(1,z)) + (1 - p_swap(:,ac_swap_index(1,z)) ) .*struct_pops_new(swap_index).per(:,ac_swap_index(1,z));
                struct_flower_one_ac.per(:,ac_swap_index(1,z)) = (1 - p_swap(:,ac_swap_index(1,z)) ).*struct_pops_new(i).per(:,ac_swap_index(1,z)) + p_swap(:,ac_swap_index(1,z)) .*struct_pops_new(swap_index).per(:,ac_swap_index(1,z));
            end
            %进行障碍物处理
            for j=1:ac_swap_num
                 %同时进行越界处理  2表示x,y坐标
                for k=1:2
                     if struct_flower_current_ac.per(k,ac_swap_index(1,j)) < pos_limit(1,1) || struct_flower_current_ac.per(k,ac_swap_index(1,j)) > pos_limit(1,2)
                        struct_flower_current_ac.per(k,ac_swap_index(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                     
                     if struct_flower_one_ac.per(k,ac_swap_index(1,j)) < pos_limit(1,1) || struct_flower_one_ac.per(k,ac_swap_index(1,j)) > pos_limit(1,2)
                        struct_flower_one_ac.per(k,ac_swap_index(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %处理障碍物  其中一个变异个体
                %左上
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point(2,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k1 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b1) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k1 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %左下
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point(3,2))
                     if (k2 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b2) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2表示相加
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k2 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b2);
                     end
                 end


                 %处理第三个区域
                 %右上角
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point(5,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point(6,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k3 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b3) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k3 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b3);
                     end
                 end

                 %处理第四区域
                 %右下角
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point(8,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point(7,2))
                     if (k4 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b4) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k4 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b4);
                     end
                 end

                %处理菱形
                %画四个区域划分
                %处理左上斜边
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,4)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k5 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b5) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%相加%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k5 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %右上
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k6 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b6) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2表示相加
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k6 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b6);
                     end
                 end


                 %处理第三个区域
                 %右下角
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,3)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,2))
                     if (k7 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b7) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k7 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b7);
                     end
                 end

                 %处理第四区域
                 %左下角
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,3)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,4))
                     if (k8 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b8) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k8 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b8);
                     end
                 end
                 
				 
                %另一个变异个体
                %处理障碍物
                %左上
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point(2,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k1 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b1) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))%左斜边
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k1 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b1);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %左下
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point(3,2))
                     if (k2 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b2) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2表示相加
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k2 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b2);
                     end
                 end


                 %处理第三个区域
                 %右上角
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point(5,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point(6,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k3 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b3) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k3 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b3);
                     end
                 end

                 %处理第四区域
                 %右下角
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point(8,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point(7,2))
                     if (k4 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b4) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k4 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b4);
                     end
                 end

                %处理菱形
                %画四个区域划分
                %处理左上斜边
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,4)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k5 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b5) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))%左斜边
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%相加%用于标记求平行线时，到底是相加还是相减的问题  1表示相减
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k5 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b5);%让其在平行的那条斜线上 
                     end
                 end


                 %处理第二个区域
                 %右上
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k6 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b6) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2表示相加
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k6 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b6);
                     end
                 end


                 %处理第三个区域
                 %右下角
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,3)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,2))
                     if (k7 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b7) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k7 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b7);
                     end
                 end

                 %处理第四区域
                 %左下角
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,3)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,4))
                     if (k8 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b8) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k8 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b8);
                     end
                 end
                
   
            end
            
            
            
            %得到新的适应值
            ac_orign_fitness = ifa_function_fitness(i,1);
            one_ac_orign_fitness = ifa_function_fitness(swap_index,1);
            
            [temp_ac_cover_fitness, temp_ac_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_flower_current_ac.per,struct_flower_current_ac.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_ac_energy_fitness] = get_energy_consume(struct_first_init.per,struct_flower_current_ac.per,struct_flower_current_ac.radius,struct_flower_current_ac.energy_init);
            fitness_temp_function_ac = weight_rate(1,1) * temp_ac_cover_fitness + weight_rate(1,2)*(1-temp_ac_energy_fitness) + weight_rate(1,3) * (1 - temp_ac_waste_fitness);

            %%判断它的适应值是否超过了之前的  如果是则进行替换
            if fitness_temp_function_ac > ac_orign_fitness
               
                struct_pops_new(i) = struct_flower_current_ac;%进行替换
                ifa_function_fitness(i,1) = fitness_temp_function_ac;
            end

            %另一个交叉的
            [temp_one_ac_cover_fitness, temp_one_ac_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_flower_one_ac.per,struct_flower_one_ac.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % 个体当前适应度
            [~,temp_one_ac_energy_fitness] = get_energy_consume(struct_first_init.per,struct_flower_one_ac.per,struct_flower_one_ac.radius,struct_flower_one_ac.energy_init);
            fitness_temp_function_one_ac = weight_rate(1,1) * temp_one_ac_cover_fitness + weight_rate(1,2)*(1-temp_one_ac_energy_fitness) + weight_rate(1,3) * (1 - temp_one_ac_waste_fitness);

            %%判断它的适应值是否超过了之前的  如果是则进行替换
            if fitness_temp_function_one_ac > one_ac_orign_fitness
                
                struct_pops_new(swap_index) = struct_flower_one_ac;%进行替换
                ifa_function_fitness(swap_index,1) = fitness_temp_function_one_ac;
            end
            
        end
    end
        

   struct_pops = struct_pops_new;%更新种群
   iter = iter+1;
end

% %写入到xlsx表中
% record_ger = 1 - record_ger;
disp('数据写入中...');
xlswrite('ifa_cover_40serson_1.xlsx',record_ger,1);


%写入到.mat文件中



figure(2);
plot(record_ger);
title('每代中最优的适应值')

figure(3);
plot(record_pop_ave);
title('每个种群的平均适应值')

figure(4);
draw_circle(struct_best_indivi.per(1,:),struct_best_indivi.per(2,:),struct_best_indivi.radius);
title('优化后的部署图');


figure(5);
draw_match(struct_first_init.per,struct_best_indivi.per,struct_first_init.radius);
title('部署好后的指派图');
hold on;





%判断连通性
[is_connec,adjacencyMatrix,adjacencyMatrix_dis ] = get_connection(struct_best_indivi.per,struct_best_indivi.radius);

if is_connec == 1
    disp('连通');
else
    disp('非连通');
end

%画最小生成树
figure(6);
draw_MST(struct_best_indivi.per(1,:),struct_best_indivi.per(2,:),struct_best_indivi.radius,adjacencyMatrix, adjacencyMatrix_dis)
title('最小生成树');
hold on;

%带回来重要消息
[match,energy,energy_rate,dis_match_first_best,energy_send_receive_match,energy_precep_match] =  get_energy_consume_end(struct_first_init.per,struct_best_indivi.per,struct_first_init.radius,struct_first_init.energy_init);

%最后匹配好后的最终位置
end_ifa_match = struct_first_init;

%匹配、还原
for i=1:N
    %位置更新
    end_ifa_match.per(:,i) = struct_best_indivi.per(:,match(i,2));
    %能量更新  距离的没算进去  发送接收的算进去了  感知的没算进去
    end_ifa_match.energy_end(1,i) = end_ifa_match.energy_init(1,i) - (0.0002*dis_match_first_best(i,2) + energy_send_receive_match(i,2) + 0.0003*energy_precep_match(i,2)); 
end
disp(['初始化时总能量：',num2str(sum(end_ifa_match.energy_init(1,:)))]);
disp(['部署后移动的总距离：',num2str(sum(dis_match_first_best(:,2)))]);
disp(['部署后移动距离消耗的能量：',num2str(0.0002*sum(dis_match_first_best(:,2)))]);
disp(['部署后消耗的总能量：',num2str(sum(end_ifa_match.energy_init(1,:)) - sum(end_ifa_match.energy_end(1,:)))]);
figure(7);
draw_circle(end_ifa_match.per(1,:),end_ifa_match.per(2,:),end_ifa_match.radius);
title('节点移动到所处的位置后部署图');
hold on;


%只为了配合获取最小生成树
[is_connec,adjacencyMatrix,adjacencyMatrix_dis ] = get_connection(end_ifa_match.per,end_ifa_match.radius);

figure(8);
draw_MST(end_ifa_match.per(1,:),end_ifa_match.per(2,:),end_ifa_match.radius,adjacencyMatrix, adjacencyMatrix_dis)
title('最后的最小生成树');
hold on;


disp(['所求函数率：',num2str(struct_best_indivi_fitness_all.function_rate)]);
disp(['对应的覆盖率：',num2str(struct_best_indivi_fitness_all.cover_rate)]);
disp(['对应的能量消耗率：',num2str(struct_best_indivi_fitness_all.energy_rate)]);
disp(['对应的节点浪费率：',num2str(struct_best_indivi_fitness_all.waste_rate)]);

%保存huiLangsousuo1的最终部署结果
end_ifa_match1 = end_ifa_match;
save end_ifa_match1.mat  end_ifa_match1;
