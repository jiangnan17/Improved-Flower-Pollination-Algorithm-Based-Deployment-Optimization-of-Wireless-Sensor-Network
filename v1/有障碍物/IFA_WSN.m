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

p=0.8;%�ж��Ƿ���ȫ���Ż����Ǿֲ��Ż�

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

struct_pops_temp =  repmat(struct_pop_per,[1 sizepop]);%��ʱ��һ����Ⱥ

energy_init_arr = zeros(1,N);
energy_end_arr = zeros(1,N);
radius_arr = zeros(1,N);
%������ε��ĸ���
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

load struct_pop_public.mat;%���ظ���Ⱥ
struct_pops = struct_pop_public;%�õ���Ⱥ����

load struct_first_init_public.mat%�����ʼ��һ����������
struct_first_init = struct_first_init_public;%�õ���ʼ����������


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
Grid_cen_x_and_y = zeros(L,W,2);%��2500���������ģ�����ÿ������������x,y
for i=1:L/1
    for j=1:W/1
        Grid_cen_x_and_y(i,j,1) = Grid_cen_x(j);%1����x
        Grid_cen_x_and_y(i,j,2) = Grid_cen_y(i);%��y����ŵ��ڶ���
    end
end


x_pos = struct_first_init.per(1,:);%��һ������   ���Ӽ��ǽ�  ��x����
y_pos = struct_first_init.per(2,:);%��һ������   ���Ӽ��ǽ�  ��y����
sersors_r = struct_first_init.radius;%��һ������ 




%���л�ͼ
figure(1);
draw_circle(x_pos,y_pos,sersors_r);
title('��ʼ������ͼ');
hold on;





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

ifa_cover_fitness = zeros(sizepop,1);%���帲����
ifa_waste_fitness = zeros(sizepop,1);%�����˷���
ifa_energy_fitness = zeros(sizepop,1);%�����ܺ���  ����1/xx��ʽ 
ifa_function_fitness = zeros(sizepop,1);%���������Ľ��
weight_rate = [1,0,0];%����Ȩ�أ�������Ӧ��������


%% Ⱥ�����
iter = 1;
record_ger = zeros(ger, 1);          % ��¼ÿ�ε�������õ���Ӧֵ 
record_pop_ave = zeros(ger,1);       % ��¼��Ⱥ��Ӧֵ��ƽ��ֵ





while iter <= ger

    %�������Ⱥ��ÿ���������Ӧֵ

    %�������Ӧֵ   �ǵ���������������  ���Եû���������
    for k=1:sizepop
        sensor_mat(1,:) = struct_pops(k).per(1,:);
        sensor_mat(2,:) = struct_pops(k).per(2,:);
        [ifa_cover_fitness(k,1), ifa_waste_fitness(k,1)] = get_Grid_cover_unit_and_rate_waste(sensor_mat,struct_pops(k).radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % ���嵱ǰ��Ӧ��
        [~,ifa_energy_fitness(k,1)] = get_energy_consume(struct_first_init.per,struct_pops(k).per,struct_pops(k).radius,struct_pops(k).energy_init);
        ifa_function_fitness(k,1) = weight_rate(1,1) * ifa_cover_fitness(k,1) + weight_rate(1,2)*(1-ifa_energy_fitness(k,1)) + weight_rate(1,3) * (1 - ifa_waste_fitness(k,1));
    end
    
   
    
    [ifa_function_fitness_sort,order_index] = sort(ifa_function_fitness);%����Ӧֵ��������
    disp(ifa_function_fitness_sort);%������ʱ�Ĵ�ӡ�Ѿ��ź���ŵ���Ӧֵ ���ڲ鿴
    

    %����������Ӧֵ�����Ÿ���
    if ifa_function_fitness(order_index(sizepop,1),1) > best_fitness
        
        %�ýṹ�屣����Ӧ�ĸ�����Ӧֵ  �Ѿ��ź����� ���������һ��
        struct_best_indivi_fitness_all.cover_rate = ifa_cover_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.waste_rate = ifa_waste_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.energy_rate = ifa_energy_fitness(order_index(sizepop,1),1);
        struct_best_indivi_fitness_all.function_rate = ifa_function_fitness(order_index(sizepop,1),1);
        
        %����������ĸ���
        struct_best_indivi = struct_pops(order_index(sizepop,1));
        %����������Ӧֵ
        best_fitness = ifa_function_fitness(order_index(sizepop,1),1);
    end
    
    record_ger(iter,1) = best_fitness;%����ô������ŵ���Ӧ��ֵ
    
    sum_fitness = sum(ifa_function_fitness);
    record_pop_ave(iter,1) = sum_fitness/sizepop;%�������Ⱥ��ƽ����Ӧ��
    
    %��һ����ʱ��Ⱥȥ����
    struct_pops_new = struct_pops;

    
%   %��������Ĵ���
    num_swap = 1;
    for j=1:sizepop
        rand_index_swap_best = randperm(N,num_swap);
        %����õĸ��������ȡһ���ֽڵ�
        for z=1:num_swap
            struct_pops_new(j).per(:,rand_index_swap_best(1,z)) = struct_best_indivi.per(:,rand_index_swap_best(1,z));
        end
    end
    
    
    
    b = 1-(1-(((ger-iter)/ger).^2)).^0.5;
    for i=1:sizepop 
        fitness_current = ifa_function_fitness(i,1);%�õ���ǰ�������Ӧֵ
         if rand(1,1) < p               %�ֲ�����
            index_rand = randperm(sizepop,2);%���ѡ������
            epsilon = b*rand(2,N);%�൱���Ǳ仯�Ĳ�����
            %%�õ��µ��м���
            struct_temp_indivi_lo = struct_pops_new(i);
            index_N_serson = randperm(N,N);%���н����д���
            for j=1:N
                struct_temp_indivi_lo.per(:,index_N_serson(1,j)) = struct_pops_new(i).per(:,index_N_serson(1,j)) + epsilon(:,j) .*(struct_pops_new(index_rand(1,1)).per(:,index_N_serson(1,j))- struct_pops_new(index_rand(1,2)).per(:,index_N_serson(1,j)));
                
                
                %ͬʱ����Խ�紦��  2��ʾx,y����
                for k=1:2
                     if struct_temp_indivi_lo.per(k,index_N_serson(1,j)) < pos_limit(1,1) || struct_temp_indivi_lo.per(k,index_N_serson(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_lo.per(k,index_N_serson(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %�����ϰ���
                %����
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%��б��
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k1 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b1);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2��ʾ���
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k2 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b2);
                     end
                 end


                 %�������������
                 %���Ͻ�
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(5,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point(6,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=50)
                     if (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k3 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b3);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point(8,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=50) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=0&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k4 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b4);
                     end
                 end

                %��������
                %���ĸ����򻮷�
                %��������б��
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,4)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,1)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,4)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,1))
                     if (k5 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b5) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))%��б��
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%���%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k5 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b5);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,1)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,2)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,1))
                     if (k6 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b6) >= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2��ʾ���
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k6 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b6);
                     end
                 end


                 %�������������
                 %���½�
                 if (struct_temp_indivi_lo.per(1,index_N_serson(1,j))>=point_diamond(1,3)&&struct_temp_indivi_lo.per(1,index_N_serson(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_lo.per(2,index_N_serson(1,j))>=point_diamond(2,3)&&struct_temp_indivi_lo.per(2,index_N_serson(1,j))<=point_diamond(2,2))
                     if (k7 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b7) <= struct_temp_indivi_lo.per(2,index_N_serson(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_temp_indivi_lo.radius(1,index_N_serson(1,j)));
                         struct_temp_indivi_lo.per(2,index_N_serson(1,j)) = (k7 * struct_temp_indivi_lo.per(1,index_N_serson(1,j)) + b7);
                     end
                 end

                 %�����������
                 %���½�
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
            
            
            %����õĸ�������´���
            lo_swap_num = 1;
            lo_serson_index = randperm(N,lo_swap_num);
            for k=1:lo_swap_num
                struct_temp_indivi_lo.per(:,lo_serson_index(1,k)) = struct_best_indivi.per(:,lo_serson_index(1,k));
            end
           
            %�ֲ������õ�����Ӧֵ

             %�õ��µ���Ӧֵ
            [temp_lo_cover_fitness, temp_lo_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_temp_indivi_lo.per,struct_temp_indivi_lo.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % ���嵱ǰ��Ӧ��
            [~,temp_lo_energy_fitness] = get_energy_consume(struct_first_init.per,struct_temp_indivi_lo.per,struct_temp_indivi_lo.radius,struct_temp_indivi_lo.energy_init);
            fitness_temp_function_lo = weight_rate(1,1) * temp_lo_cover_fitness + weight_rate(1,2)*(1-temp_lo_energy_fitness) + weight_rate(1,3) * (1 - temp_lo_waste_fitness);

            %%�ж�������Ӧֵ�Ƿ񳬹���֮ǰ��  �����������滻
            if fitness_temp_function_lo > fitness_current
                struct_pops_new(i) = struct_temp_indivi_lo;%�����滻
                
            end
         else%����ȫ�ֵ�����
            L1 = Levy2(dimension,iter,ger);%�õ�����
            index_rand = randperm(N,N);%���ѡ������
            %%�õ��µ��м���
            struct_temp_indivi_glo = struct_pops_new(i);
            for j=1:N
                struct_temp_indivi_glo.per(:,index_rand(1,j)) = struct_pops_new(i).per(:,index_rand(1,j)) + L1' .*(struct_pops_new(i).per(:,index_rand(1,j))- struct_best_indivi.per(:,index_rand(1,j)));%�õ���ʱ����
            
                 %ͬʱ����Խ�紦��  2��ʾx,y����
                for k=1:2
                     if struct_temp_indivi_glo.per(k,index_rand(1,j)) < pos_limit(1,1) || struct_temp_indivi_glo.per(k,index_rand(1,j)) > pos_limit(1,2)
                        struct_temp_indivi_glo.per(k,index_rand(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %�����ϰ���
                %����
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1) <= struct_temp_indivi_glo.per(2,index_rand(1,j))%��б��
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k1 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b1);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(3,2))
                     if (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2��ʾ���
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k2 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b2);
                     end
                 end


                 %�������������
                 %���Ͻ�
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(5,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point(6,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=50)
                     if (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k3 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b3);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point(8,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=50) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=0&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point(7,2))
                     if (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k4 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b4);
                     end
                 end

                %��������
                %���ĸ����򻮷�
                %��������б��
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,4)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,1)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,4)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,1))
                     if (k5 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b5) >= struct_temp_indivi_glo.per(2,index_rand(1,j))%��б��
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%���%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k5 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b5);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,1)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,2)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,1))
                     if (k6 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b6) >= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2��ʾ���
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k6 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b6);
                     end
                 end


                 %�������������
                 %���½�
                 if (struct_temp_indivi_glo.per(1,index_rand(1,j))>=point_diamond(1,3)&&struct_temp_indivi_glo.per(1,index_rand(1,j))<=point_diamond(1,2)) && (struct_temp_indivi_glo.per(2,index_rand(1,j))>=point_diamond(2,3)&&struct_temp_indivi_glo.per(2,index_rand(1,j))<=point_diamond(2,2))
                     if (k7 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b7) <= struct_temp_indivi_glo.per(2,index_rand(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_temp_indivi_glo.radius(1,index_rand(1,j)));
                         struct_temp_indivi_glo.per(2,index_rand(1,j)) = (k7 * struct_temp_indivi_glo.per(1,index_rand(1,j)) + b7);
                     end
                 end

                 %�����������
                 %���½�
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
            
            %����õĸ�������´���
            glo_swap_num = 2;
            glo_serson_index = randperm(N,glo_swap_num);
            for k=1:glo_swap_num
                struct_temp_indivi_glo.per(:,glo_serson_index(1,k)) = struct_best_indivi.per(:,glo_serson_index(1,k));
            end
            
             %�õ��µ���Ӧֵ
            [temp_glo_cover_fitness, temp_glo_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_temp_indivi_glo.per,struct_temp_indivi_glo.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % ���嵱ǰ��Ӧ��
            [~,temp_glo_energy_fitness] = get_energy_consume(struct_first_init.per,struct_temp_indivi_glo.per,struct_temp_indivi_glo.radius,struct_temp_indivi_glo.energy_init);
            fitness_temp_function_glo = weight_rate(1,1) * temp_glo_cover_fitness + weight_rate(1,2)*(1-temp_glo_energy_fitness) + weight_rate(1,3) * (1 - temp_glo_waste_fitness);

            %%�ж�������Ӧֵ�Ƿ񳬹���֮ǰ��  �����������滻
            if fitness_temp_function_glo > fitness_current
                struct_pops_new(i) = struct_temp_indivi_glo;%�����滻
                
            end
           
         end
        
        %�ֲ��Ľ���
        if rand(1,1) < 0.4;
            rand_swap_p = b*rand(1,1);%���ڽ�����Ϣ
            flower_rand_one_va = randperm(sizepop,2);
            %�õ�һ��������flower���н���
            for kk=1:2
                if i ~= flower_rand_one_va(1,kk)
                    swap_index = flower_rand_one_va(1,kk);
                    break;
                end
            end
            %���н���
            %����
            p_swap = b*rand(2,N);
            struct_flower_current_ac = struct_pops_new(i);
            struct_flower_one_ac = struct_pops_new(swap_index);
            ac_swap_num = 1;
            ac_swap_index = randperm(N,ac_swap_num);
            %����
            for z=1:ac_swap_num
               
                struct_flower_current_ac.per(:,ac_swap_index(1,z)) = p_swap(:,ac_swap_index(1,z)) .* struct_pops_new(i).per(:,ac_swap_index(1,z)) + (1 - p_swap(:,ac_swap_index(1,z)) ) .*struct_pops_new(swap_index).per(:,ac_swap_index(1,z));
                struct_flower_one_ac.per(:,ac_swap_index(1,z)) = (1 - p_swap(:,ac_swap_index(1,z)) ).*struct_pops_new(i).per(:,ac_swap_index(1,z)) + p_swap(:,ac_swap_index(1,z)) .*struct_pops_new(swap_index).per(:,ac_swap_index(1,z));
            end
            %�����ϰ��ﴦ��
            for j=1:ac_swap_num
                 %ͬʱ����Խ�紦��  2��ʾx,y����
                for k=1:2
                     if struct_flower_current_ac.per(k,ac_swap_index(1,j)) < pos_limit(1,1) || struct_flower_current_ac.per(k,ac_swap_index(1,j)) > pos_limit(1,2)
                        struct_flower_current_ac.per(k,ac_swap_index(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                     
                     if struct_flower_one_ac.per(k,ac_swap_index(1,j)) < pos_limit(1,1) || struct_flower_one_ac.per(k,ac_swap_index(1,j)) > pos_limit(1,2)
                        struct_flower_one_ac.per(k,ac_swap_index(1,j)) = pos_limit(1,1) + rand(1,1)*(pos_limit(1,2) - pos_limit(1,1));
                     end
                end
                
                %�����ϰ���  ����һ���������
                %����
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point(2,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k1 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b1) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))%��б��
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k1 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b1);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point(3,2))
                     if (k2 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b2) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2��ʾ���
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k2 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b2);
                     end
                 end


                 %�������������
                 %���Ͻ�
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point(5,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point(6,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k3 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b3) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k3 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b3);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point(8,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point(7,2))
                     if (k4 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b4) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k4 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b4);
                     end
                 end

                %��������
                %���ĸ����򻮷�
                %��������б��
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,1)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,4)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k5 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b5) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))%��б��
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%���%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k5 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b5);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,1)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,2)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k6 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b6) >= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2��ʾ���
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k6 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b6);
                     end
                 end


                 %�������������
                 %���½�
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,3)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,2))
                     if (k7 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b7) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k7 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b7);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_flower_current_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_current_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,3)) && (struct_flower_current_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_current_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,4))
                     if (k8 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b8) <= struct_flower_current_ac.per(2,ac_swap_index(1,j))
                         k8 = 1;
                         b8 = 40;
                         case_b = 1;
                         [k8,b8] = get_new_function(k8,b8,case_b,struct_flower_current_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_current_ac.per(2,ac_swap_index(1,j)) = (k8 * struct_flower_current_ac.per(1,ac_swap_index(1,j)) + b8);
                     end
                 end
                 
				 
                %��һ���������
                %�����ϰ���
                %����
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point(2,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k1 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b1) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))%��б��
                         k1 = 1;
                         b1 = 35;
                         case_b = 1;%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k1,b1] = get_new_function(k1,b1,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k1 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b1);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point(3,2))
                     if (k2 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b2) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k2 = -1;
                         b2 = 15;
                         case_b = 2;%2��ʾ���
                         [k2,b2] = get_new_function(k2,b2,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k2 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b2);
                     end
                 end


                 %�������������
                 %���Ͻ�
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point(5,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point(6,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=50)
                     if (k3 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b3) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k3 = -1;
                         b3 = 85;
                         case_b = 1;
                         [k3,b3] = get_new_function(k3,b3,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k3 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b3);
                     end
                 end

                 %�����������
                 %���½�
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point(8,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=50) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=0&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point(7,2))
                     if (k4 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b4) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k4 = 1;
                         b4 = -35;
                         case_b = 2;
                         [k4,b4] = get_new_function(k4,b4,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k4 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b4);
                     end
                 end

                %��������
                %���ĸ����򻮷�
                %��������б��
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,4)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,1)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,4)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k5 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b5) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))%��б��
                         k5 = 1;
                         b5 = 10;
                         case_b = 2;%���%���ڱ����ƽ����ʱ����������ӻ������������  1��ʾ���
                         [k5,b5] = get_new_function(k5,b5,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k5 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b5);%������ƽ�е�����б���� 
                     end
                 end


                 %����ڶ�������
                 %����
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,1)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,2)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,1))
                     if (k6 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b6) >= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k6 = -1;
                         b6 = 60;
                         case_b = 2;%2��ʾ���
                         [k6,b6] = get_new_function(k6,b6,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k6 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b6);
                     end
                 end


                 %�������������
                 %���½�
                 if (struct_flower_one_ac.per(1,ac_swap_index(1,j))>=point_diamond(1,3)&&struct_flower_one_ac.per(1,ac_swap_index(1,j))<=point_diamond(1,2)) && (struct_flower_one_ac.per(2,ac_swap_index(1,j))>=point_diamond(2,3)&&struct_flower_one_ac.per(2,ac_swap_index(1,j))<=point_diamond(2,2))
                     if (k7 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b7) <= struct_flower_one_ac.per(2,ac_swap_index(1,j))
                         k7 = 1;
                         b7 = -10;
                         case_b = 1;
                         [k7,b7] = get_new_function(k7,b7,case_b,struct_flower_one_ac.radius(1,ac_swap_index(1,j)));
                         struct_flower_one_ac.per(2,ac_swap_index(1,j)) = (k7 * struct_flower_one_ac.per(1,ac_swap_index(1,j)) + b7);
                     end
                 end

                 %�����������
                 %���½�
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
            
            
            
            %�õ��µ���Ӧֵ
            ac_orign_fitness = ifa_function_fitness(i,1);
            one_ac_orign_fitness = ifa_function_fitness(swap_index,1);
            
            [temp_ac_cover_fitness, temp_ac_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_flower_current_ac.per,struct_flower_current_ac.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % ���嵱ǰ��Ӧ��
            [~,temp_ac_energy_fitness] = get_energy_consume(struct_first_init.per,struct_flower_current_ac.per,struct_flower_current_ac.radius,struct_flower_current_ac.energy_init);
            fitness_temp_function_ac = weight_rate(1,1) * temp_ac_cover_fitness + weight_rate(1,2)*(1-temp_ac_energy_fitness) + weight_rate(1,3) * (1 - temp_ac_waste_fitness);

            %%�ж�������Ӧֵ�Ƿ񳬹���֮ǰ��  �����������滻
            if fitness_temp_function_ac > ac_orign_fitness
               
                struct_pops_new(i) = struct_flower_current_ac;%�����滻
                ifa_function_fitness(i,1) = fitness_temp_function_ac;
            end

            %��һ�������
            [temp_one_ac_cover_fitness, temp_one_ac_waste_fitness] = get_Grid_cover_unit_and_rate_waste(struct_flower_one_ac.per,struct_flower_one_ac.radius(1,:),per_sersons_num(1,:),per_sersons_radius_type(1,:)) ; % ���嵱ǰ��Ӧ��
            [~,temp_one_ac_energy_fitness] = get_energy_consume(struct_first_init.per,struct_flower_one_ac.per,struct_flower_one_ac.radius,struct_flower_one_ac.energy_init);
            fitness_temp_function_one_ac = weight_rate(1,1) * temp_one_ac_cover_fitness + weight_rate(1,2)*(1-temp_one_ac_energy_fitness) + weight_rate(1,3) * (1 - temp_one_ac_waste_fitness);

            %%�ж�������Ӧֵ�Ƿ񳬹���֮ǰ��  �����������滻
            if fitness_temp_function_one_ac > one_ac_orign_fitness
                
                struct_pops_new(swap_index) = struct_flower_one_ac;%�����滻
                ifa_function_fitness(swap_index,1) = fitness_temp_function_one_ac;
            end
            
        end
    end
        

   struct_pops = struct_pops_new;%������Ⱥ
   iter = iter+1;
end

% %д�뵽xlsx����
% record_ger = 1 - record_ger;
disp('����д����...');
xlswrite('ifa_cover_40serson_1.xlsx',record_ger,1);


%д�뵽.mat�ļ���



figure(2);
plot(record_ger);
title('ÿ�������ŵ���Ӧֵ')

figure(3);
plot(record_pop_ave);
title('ÿ����Ⱥ��ƽ����Ӧֵ')

figure(4);
draw_circle(struct_best_indivi.per(1,:),struct_best_indivi.per(2,:),struct_best_indivi.radius);
title('�Ż���Ĳ���ͼ');


figure(5);
draw_match(struct_first_init.per,struct_best_indivi.per,struct_first_init.radius);
title('����ú��ָ��ͼ');
hold on;





%�ж���ͨ��
[is_connec,adjacencyMatrix,adjacencyMatrix_dis ] = get_connection(struct_best_indivi.per,struct_best_indivi.radius);

if is_connec == 1
    disp('��ͨ');
else
    disp('����ͨ');
end

%����С������
figure(6);
draw_MST(struct_best_indivi.per(1,:),struct_best_indivi.per(2,:),struct_best_indivi.radius,adjacencyMatrix, adjacencyMatrix_dis)
title('��С������');
hold on;

%��������Ҫ��Ϣ
[match,energy,energy_rate,dis_match_first_best,energy_send_receive_match,energy_precep_match] =  get_energy_consume_end(struct_first_init.per,struct_best_indivi.per,struct_first_init.radius,struct_first_init.energy_init);

%���ƥ��ú������λ��
end_ifa_match = struct_first_init;

%ƥ�䡢��ԭ
for i=1:N
    %λ�ø���
    end_ifa_match.per(:,i) = struct_best_indivi.per(:,match(i,2));
    %��������  �����û���ȥ  ���ͽ��յ����ȥ��  ��֪��û���ȥ
    end_ifa_match.energy_end(1,i) = end_ifa_match.energy_init(1,i) - (0.0002*dis_match_first_best(i,2) + energy_send_receive_match(i,2) + 0.0003*energy_precep_match(i,2)); 
end
disp(['��ʼ��ʱ��������',num2str(sum(end_ifa_match.energy_init(1,:)))]);
disp(['������ƶ����ܾ��룺',num2str(sum(dis_match_first_best(:,2)))]);
disp(['������ƶ��������ĵ�������',num2str(0.0002*sum(dis_match_first_best(:,2)))]);
disp(['��������ĵ���������',num2str(sum(end_ifa_match.energy_init(1,:)) - sum(end_ifa_match.energy_end(1,:)))]);
figure(7);
draw_circle(end_ifa_match.per(1,:),end_ifa_match.per(2,:),end_ifa_match.radius);
title('�ڵ��ƶ���������λ�ú���ͼ');
hold on;


%ֻΪ����ϻ�ȡ��С������
[is_connec,adjacencyMatrix,adjacencyMatrix_dis ] = get_connection(end_ifa_match.per,end_ifa_match.radius);

figure(8);
draw_MST(end_ifa_match.per(1,:),end_ifa_match.per(2,:),end_ifa_match.radius,adjacencyMatrix, adjacencyMatrix_dis)
title('������С������');
hold on;


disp(['�������ʣ�',num2str(struct_best_indivi_fitness_all.function_rate)]);
disp(['��Ӧ�ĸ����ʣ�',num2str(struct_best_indivi_fitness_all.cover_rate)]);
disp(['��Ӧ�����������ʣ�',num2str(struct_best_indivi_fitness_all.energy_rate)]);
disp(['��Ӧ�Ľڵ��˷��ʣ�',num2str(struct_best_indivi_fitness_all.waste_rate)]);

%����huiLangsousuo1�����ղ�����
end_ifa_match1 = end_ifa_match;
save end_ifa_match1.mat  end_ifa_match1;
