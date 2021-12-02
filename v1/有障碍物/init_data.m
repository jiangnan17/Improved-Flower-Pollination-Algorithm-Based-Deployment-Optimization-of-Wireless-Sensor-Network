%%������
clc;
clear ;
close all;
%ɾ����Ӧ���ļ�


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
N = 25;%30���������ڵ�
sizepop = 50;%��Ⱥ��ģ
dimension = 2;% �ռ�ά��  ǰ�з�x��y�������зŰ뾶
ger = 10;% ����������
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

%��������ε��ĸ���
syms x y;%�ȶ���һ������
%���Ͻ�
k1 = 1;
b1 = 35;
x1_up = solve(k1*x+b1==50,x);%���Ͻǵ�б�ߵ��ϸ�����
y1_down = solve(k1*0+b1==y,y);

%���½�
k2 = -1;
b2 = 15;
y2_up = solve(k2*0+b2==y,y);
x2_down = solve(k2*x+b2==0,x);


%���Ͻ�
k3 = -1;
b3 = 85;
x3_up = solve(k3*x+b3==50,x);
y3_down = solve(k3*50+b3==y,y);

%���½�
k4 = 1;
b4 = -35;
y4_up = solve(k4*50+b4==y,y);
x4_down = solve(k4*x+b4==0,x);

%����������֤��ϣ���ȫ��ȷ
point = zeros(8,2);%�洢��Щ��  ����  ��������
point(1,:) = [x1_up,50];
point(2,:) = [0,y1_down];
point(3,:) = [0,y2_up];
point(4,:) = [x2_down,0];
point(5,:) = [x3_up,50];
point(6,:) = [50,y3_down];
point(7,:) = [50,y4_up];
point(8,:) = [x4_down,0];

%���εļ���
point_diamond = zeros(2,4);%���ε��ĸ��㣬��λ��˳ʱ�� ��һ��Ϊ�� ����Ϊ��
%����������ε��ĸ���
syms x y;%�ȶ���һ������
%���Ͻ�
k5 = 1;
b5 = 10;
%���ʲô������  ֱ�ӿ��Կ�����

point_diamond(1,1) = 25;
point_diamond(2,1) = 35;

%���Ͻ�
k6 = -1;
b6 = 60;

point_diamond(1,2) = 35;
point_diamond(2,2) = 25;

%���½�
k7 = 1;
b7 = -10;

point_diamond(1,3) = 25;
point_diamond(2,3) = 15;

%���½�
k8 = 1;
b8 = 40;

point_diamond(1,4) = 15;
point_diamond(2,4) = 25;

%���ɽڵ����ݡ��������뾶�ṹ��  ͬʱ����������䵽�ϰ�������
for i=1:sizepop
    
    %��ʼ���뾶������  ��Ϊ�����漰���뾶��ʹ��  �����ȳ�ʼ���뾶
    struct_pops(i).radius(1,:) = radius_arr;
    struct_pops(i).energy_init(1,:) = energy_init_arr;
    struct_pops(i).energy_end(1,:) = energy_end_arr;
    
    %��ʼ���ڵ����
    struct_pops(i).sersons_num = [r_max_num,r_mid_num,r_min_num];
    
    %��ʼ���ڵ�λ������
    x = L*rand(1,N);%������ɽڵ�ĺ�����
    y = W*rand(1,N);%������ɽڵ��������
    struct_pops(i).per(1,:) = x;
    struct_pops(i).per(2,:) = y;
    for j=1:N
        %���ĸ����򻮷�
        %��������б��
         if (struct_pops(i).per(1,j)>=0&&struct_pops(i).per(1,j)<=point(1,1)) && (struct_pops(i).per(2,j)>=point(2,2)&&struct_pops(i).per(2,j)<=50)
             if (k1 * struct_pops(i).per(1,j) + b1) <= struct_pops(i).per(2,j)%��б��
                 k1 = 1;
                 b1 = 35;
                 case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                 [k1,b1] = get_new_function(k1,b1,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k1 * struct_pops(i).per(1,j) + b1);%������ƽ�е�����б���� 
             end
         end
         
         
         %����ڶ�������
         %����
         if (struct_pops(i).per(1,j)>=0&&struct_pops(i).per(1,j)<=point(1,1)) && (struct_pops(i).per(2,j)>=0&&struct_pops(i).per(2,j)<=point(3,2))
             if (k2 * struct_pops(i).per(1,j) + b2) >= struct_pops(i).per(2,j)
                 k2 = -1;
                 b2 = 15;
                 case_b = 2;%2��ʾ���
                 [k2,b2] = get_new_function(k2,b2,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k2 * struct_pops(i).per(1,j) + b2);
             end
         end
         
         
         %�������������
         %���Ͻ�
         if (struct_pops(i).per(1,j)>=point(5,1)&&struct_pops(i).per(1,j)<=50) && (struct_pops(i).per(2,j)>=point(6,2)&&struct_pops(i).per(2,j)<=50)
             if (k3 * struct_pops(i).per(1,j) + b3) <= struct_pops(i).per(2,j)
                 k3 = -1;
                 b3 = 85;
                 case_b = 1;
                 [k3,b3] = get_new_function(k3,b3,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k3 * struct_pops(i).per(1,j) + b3);
             end
         end
         
         %�����������
         %���½�
         if (struct_pops(i).per(1,j)>=point(8,1)&&struct_pops(i).per(1,j)<=50) && (struct_pops(i).per(2,j)>=0&&struct_pops(i).per(2,j)<=point(7,2))
             if (k4 * struct_pops(i).per(1,j) + b4) >= struct_pops(i).per(2,j)
                 k4 = 1;
                 b4 = -35;
                 case_b = 2;
                 [k4,b4] = get_new_function(k4,b4,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k4 * struct_pops(i).per(1,j) + b4);
             end
         end
         
        %��������
        %���ĸ����򻮷�
        %��������б��
         if (struct_pops(i).per(1,j)>=point_diamond(1,4)&&struct_pops(i).per(1,j)<=point_diamond(1,1)) && (struct_pops(i).per(2,j)>=point_diamond(2,4)&&struct_pops(i).per(2,j)<=point_diamond(2,1))
             if (k5 * struct_pops(i).per(1,j) + b5) >= struct_pops(i).per(2,j)%��б��
                 k5 = 1;
                 b5 = 10;
                 case_b = 2;%���%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                 [k5,b5] = get_new_function(k5,b5,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k5 * struct_pops(i).per(1,j) + b5);%������ƽ�е�����б���� 
             end
         end
         
         
         %����ڶ�������
         %����
         if (struct_pops(i).per(1,j)>=point_diamond(1,1)&&struct_pops(i).per(1,j)<=point_diamond(1,2)) && (struct_pops(i).per(2,j)>=point_diamond(2,2)&&struct_pops(i).per(2,j)<=point_diamond(2,1))
             if (k6 * struct_pops(i).per(1,j) + b6) >= struct_pops(i).per(2,j)
                 k6 = -1;
                 b6 = 60;
                 case_b = 2;%2��ʾ���
                 [k6,b6] = get_new_function(k6,b6,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k6 * struct_pops(i).per(1,j) + b6);
             end
         end
         
         
         %�������������
         %���½�
         if (struct_pops(i).per(1,j)>=point_diamond(1,3)&&struct_pops(i).per(1,j)<=point_diamond(1,2)) && (struct_pops(i).per(2,j)>=point_diamond(2,3)&&struct_pops(i).per(2,j)<=point_diamond(2,2))
             if (k7 * struct_pops(i).per(1,j) + b7) <= struct_pops(i).per(2,j)
                 k7 = 1;
                 b7 = -10;
                 case_b = 1;
                 [k7,b7] = get_new_function(k7,b7,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k7 * struct_pops(i).per(1,j) + b7);
             end
         end
         
         %�����������
         %���½�
         if (struct_pops(i).per(1,j)>=point_diamond(1,4)&&struct_pops(i).per(1,j)<=point_diamond(1,3)) && (struct_pops(i).per(2,j)>=point_diamond(2,3)&&struct_pops(i).per(2,j)<=point_diamond(2,4))
             if (k8 * struct_pops(i).per(1,j) + b8) <= struct_pops(i).per(2,j)
                 k8 = 1;
                 b8 = 40;
                 case_b = 1;
                 [k8,b8] = get_new_function(k8,b8,case_b,struct_pops(i).radius(1,j));
                 struct_pops(i).per(2,j) = (k8 * struct_pops(i).per(1,j) + b8);
             end
         end
         
    end
    
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
%������������
%������������һ��Ϊ��һ��
Grid_cen_x_and_y = zeros(L,W,2);%��2500���������ģ�����ÿ������������x,y
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);%1����x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);%��y����ŵ��ڶ���
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
%Ϊ�������������µ�һ�������sizepop����������ƴ�ն���
new_first_per = zeros(2,N);
piece_together_index = randperm(N,N);
rand_indivi = randperm(sizepop,N);%ע��NС��sizepop
for i=1:N
    new_first_per(:,i) = struct_pops(rand_indivi(1,i)).per(:,piece_together_index(1,i));
end
struct_pops(1).per = new_first_per;
struct_pop_public = struct_pops;%�����ĳ�ʼ��Ⱥ
%������Ⱥ��wolf_pop_public.mat�ļ���
save struct_pop_public.mat struct_pop_public;

disp('��ʼ���������');
