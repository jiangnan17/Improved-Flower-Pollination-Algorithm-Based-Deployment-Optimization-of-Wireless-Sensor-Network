%%主程序
clc;
clear ;
close all;



global N;
global M;
global L;
global W;
global Grid_cen_x;
global Grid_cen_y;
global Grid_cen_x_and_y;
global ger;
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
N = 35;%30个传感器节点
sizepop = 50;%种群规模
dimension = 2;% 空间维数  前行放x、y，第三行放半径
ger = 200;% 最大迭代次数
pos_limit = [0, 50];            % 设置位置参数限制
%个数限制
r_max_num = 1;%序号为1-5
r_mid_num = 2;%序号为6-10
r_min_num = N - r_max_num - r_mid_num; %序号为11-N


struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);%结构体类型
struct_pops = repmat(struct_pop_per,[1 sizepop]);%生成结构体数组,1行20列

struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);%临时的一个种群

energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);

%先初始化了能量和半径  这个是统一的
for k=1:N
    if k>=1 && k <=r_max_num 
       radius_arr(1,k) = r_max; 
       energy_init_arr(1,k) = energy_max;
       energy_end_arr(1,k) = energy_max;
    elseif k>=r_max_num+1 && k <= (r_max_num + r_mid_num)
       radius_arr(1,k) = r_mid; 
       energy_init_arr(1,k) = energy_mid;
       energy_end_arr(1,k) = energy_mid;
    else
       radius_arr(1,k) = r_min; 
       energy_init_arr(1,k) = energy_min;
       energy_end_arr(1,k) = energy_min;
    end
end


%生成节点数据、能量、半径结构体
for i=1:sizepop
    %初始化节点位置坐标
    x = L*rand(1,N);%随机生成节点的横坐标
    y = W*rand(1,N);%随机生成节点的列坐标
    struct_pops(i).per(1,:) = x;
    struct_pops(i).per(2,:) = y;
    
    %初始化半径和能量
    struct_pops(i).radius(1,:) = radius_arr;
    struct_pops(i).energy_init(1,:) = energy_init_arr;
    struct_pops(i).energy_end(1,:) = energy_end_arr;
    
    %初始化节点个数
    struct_pops(i).sersons_num = [r_max_num,r_mid_num,r_min_num];
end



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
Grid_cen_x_and_y = zeros(2,M);
k = 1;
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(1,k) = Grid_cen_x(j);%把x坐标放到第一行
        Grid_cen_x_and_y(2,k) = Grid_cen_y(i);%把y坐标放到第二行
        k = k+1;%下一个点
    end
end


x_pos = struct_pops(1).per(1,:);%第一个粒子   粒子即是解  的x坐标
y_pos = struct_pops(1).per(2,:);%第一个粒子   粒子即是解  的y坐标
sersors_r = struct_pops(1).radius;%第一个粒子 




%进行画图
figure(1);
draw_circle(x_pos,y_pos,sersors_r);
title('初始化部署图');
hold on;
struct_first_init_wolf = struct_pops(1);%存第一只狼数据

%保存一份第一只狼的数据
struct_first_init_public = struct_first_init_wolf;
save struct_first_init_public.mat struct_first_init_public;



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
wolf_one_dis =  zeros(2,N);%第一只狼距离食物的距离
wolf_two_dis = zeros(2,N);%第二只狼距离食物
wolf_three_dis = zeros(2,N);%第三只狼距离食物
wolf_one = zeros(2,N);%第一只狼
wlof_two = zeros(2,N);%第二只狼
wolf_three = zeros(2,N);%第三只狼
wolf_cover_fitness = zeros(sizepop,1);%个体覆盖率
wolf_waste_fitness = zeros(sizepop,1);%个体浪费率
wolf_energy_fitness = zeros(sizepop,1);%个体能耗率  化成1/xx形式 
wolf_function_fitness = zeros(sizepop,1);%三个函数的结合
weight_rate = [1,0,0];%设置权重，覆盖率应该是最大的


wolf_fitness_temp = zeros(sizepop,1);%临时保存当前适应度
%% 群体更新
iter = 1;
record_ger = zeros(ger, 1);          % 记录每次迭代中最好的适应值 
record_pop_ave = zeros(ger,1);       % 记录种群适应值的平均值


%因为第一代中的第一个个体用来当做初始，故得考虑重新生成第一个个体 只生成坐标
new_first_per = [L*rand(1,N);W*rand(1,N)];
struct_pops(1).per = new_first_per;

struct_pop_public = struct_pops;%公共的初始种群
%保存种群到wolf_pop_public.mat文件中
save struct_pop_public.mat struct_pop_public;


disp('数据初始化完成');