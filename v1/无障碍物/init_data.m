%%������
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
L = 50;%��
W = 50;%��
%����1ƽ����һ������
M = 2500;%��������
r_max = 7;%��֪�뾶Ϊ5
r_mid = 6;
r_min = 5;
energy_max = 100;%��������
energy_mid = 90;
energy_min = 80;

per_sersons_radius_type = [r_max,r_mid,r_min];
%�������Ϊ5��ʣ��ΪС
N = 35;%30���������ڵ�
sizepop = 50;%��Ⱥ��ģ
dimension = 2;% �ռ�ά��  ǰ�з�x��y�������зŰ뾶
ger = 200;% ����������
pos_limit = [0, 50];            % ����λ�ò�������
%��������
r_max_num = 1;%���Ϊ1-5
r_mid_num = 2;%���Ϊ6-10
r_min_num = N - r_max_num - r_mid_num; %���Ϊ11-N


struct_pop_per = struct('per',[],'radius',[],'energy_init',[],'energy_end',[],'sersons_num',[]);%�ṹ������
struct_pops = repmat(struct_pop_per,[1 sizepop]);%���ɽṹ������,1��20��

struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);%��ʱ��һ����Ⱥ

energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);

%�ȳ�ʼ���������Ͱ뾶  �����ͳһ��
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


%���ɽڵ����ݡ��������뾶�ṹ��
for i=1:sizepop
    %��ʼ���ڵ�λ������
    x = L*rand(1,N);%������ɽڵ�ĺ�����
    y = W*rand(1,N);%������ɽڵ��������
    struct_pops(i).per(1,:) = x;
    struct_pops(i).per(2,:) = y;
    
    %��ʼ���뾶������
    struct_pops(i).radius(1,:) = radius_arr;
    struct_pops(i).energy_init(1,:) = energy_init_arr;
    struct_pops(i).energy_end(1,:) = energy_end_arr;
    
    %��ʼ���ڵ����
    struct_pops(i).sersons_num = [r_max_num,r_mid_num,r_min_num];
end



%%��ʼ�Ĳ����ͼ  �õ�һ��������ȥ��ʼ��ͼ

%��������������
X_mat = (0:1:50);%x����
Y_mat = (0:1:50);%y����
Grid_cen_x = zeros(1,L/1);%�������ĵ�x����
Grid_cen_y = zeros(1,W/1);%�������ĵ�y����
%ǰ���������֮�ͳ���2
for i=1:L/1
    Grid_cen_x(i) = (X_mat(i)+X_mat(i+1))/2;
    Grid_cen_y(i) = (Y_mat(i)+Y_mat(i+1))/2;
end

%%�Ѻ������궪��һ����ά������
%����ת����  ��һ�з�x�� �ڶ��з�y�ᣬͬһ�з�һ��������
%���ȴ濿��x����ĵ�һ�У�Ȼ�����ϴ�ڶ���
%������������
Grid_cen_x_and_y = zeros(2,M);
k = 1;
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(1,k) = Grid_cen_x(j);%��x����ŵ���һ��
        Grid_cen_x_and_y(2,k) = Grid_cen_y(i);%��y����ŵ��ڶ���
        k = k+1;%��һ����
    end
end


x_pos = struct_pops(1).per(1,:);%��һ������   ���Ӽ��ǽ�  ��x����
y_pos = struct_pops(1).per(2,:);%��һ������   ���Ӽ��ǽ�  ��y����
sersors_r = struct_pops(1).radius;%��һ������ 




%���л�ͼ
figure(1);
draw_circle(x_pos,y_pos,sersors_r);
title('��ʼ������ͼ');
hold on;
struct_first_init_wolf = struct_pops(1);%���һֻ������

%����һ�ݵ�һֻ�ǵ�����
struct_first_init_public = struct_first_init_wolf;
save struct_first_init_public.mat struct_first_init_public;



%���������ĵ������ŵ�������
sensor_mat = zeros(2,N);%Ԥ�����ڴ�
for i=1:N
    sensor_mat(1,i) = x_pos(i);
    sensor_mat(2,i) = y_pos(i);
end

%������ֽڵ����Ŀ
per_sersons_num = struct_pops(1).sersons_num;
%������ֽڵ�İ뾶


%%��ʼ���õ����ϸ��ʺͽڵ��˷���
[cover_rate,waste_rate] =  get_Grid_cover_unit_and_rate_waste(sensor_mat,sersors_r,per_sersons_num,per_sersons_radius_type);
disp(['��ʼ���ĸ����ʣ�',num2str(cover_rate)]);
disp(['��ʼ�����˷��ʣ�',num2str(waste_rate)]);



%%������ͨ��  ��һ���ڵ��ʼ��ʱ����ͨ��
is_connec = get_connection(sensor_mat,sersors_r);
if is_connec==1
    disp('��ͨ');
else
    disp('����ͨ');
end




%%��ʼ����Ⱥ��ʷֵΪ  ����С   infΪ�����
best_fitness = -inf;                         % ��Ⱥ��ʷ�����Ӧ��  
struct_best_indivi = struct_pop_per;                 % �����������
struct_best_indivi_fitness = struct('cover_rate',[],'waste_rate',[],'energy_rate',[],'function_rate',[]);%��Ӧֵ�ṹ������
struct_best_indivi_fitness_all = repmat(struct_best_indivi_fitness,[1 1]);%Ԥ�����ڴ�
%%�������Ѿ���ͼ��
%%����Ϊ������ʼ��ʱ��50�����ӿ�ʼ��λ��
wolf_one_dis =  zeros(2,N);%��һֻ�Ǿ���ʳ��ľ���
wolf_two_dis = zeros(2,N);%�ڶ�ֻ�Ǿ���ʳ��
wolf_three_dis = zeros(2,N);%����ֻ�Ǿ���ʳ��
wolf_one = zeros(2,N);%��һֻ��
wlof_two = zeros(2,N);%�ڶ�ֻ��
wolf_three = zeros(2,N);%����ֻ��
wolf_cover_fitness = zeros(sizepop,1);%���帲����
wolf_waste_fitness = zeros(sizepop,1);%�����˷���
wolf_energy_fitness = zeros(sizepop,1);%�����ܺ���  ����1/xx��ʽ 
wolf_function_fitness = zeros(sizepop,1);%���������Ľ��
weight_rate = [1,0,0];%����Ȩ�أ�������Ӧ��������


wolf_fitness_temp = zeros(sizepop,1);%��ʱ���浱ǰ��Ӧ��
%% Ⱥ�����
iter = 1;
record_ger = zeros(ger, 1);          % ��¼ÿ�ε�������õ���Ӧֵ 
record_pop_ave = zeros(ger,1);       % ��¼��Ⱥ��Ӧֵ��ƽ��ֵ


%��Ϊ��һ���еĵ�һ����������������ʼ���ʵÿ����������ɵ�һ������ ֻ��������
new_first_per = [L*rand(1,N);W*rand(1,N)];
struct_pops(1).per = new_first_per;

struct_pop_public = struct_pops;%�����ĳ�ʼ��Ⱥ
%������Ⱥ��wolf_pop_public.mat�ļ���
save struct_pop_public.mat struct_pop_public;


disp('���ݳ�ʼ�����');