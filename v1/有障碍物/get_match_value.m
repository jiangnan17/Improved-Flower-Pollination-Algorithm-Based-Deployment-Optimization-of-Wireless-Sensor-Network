function [match,value,matrix_dis_first_any] = get_match_value(first_init_wolf,any_indivi,sensor_r)

    [~,N] = size(first_init_wolf);%�õ��ж��ٸ��ڵ�
    
    num_7 = 0;%�뾶Ϊ7�ĸ���
    num_6 = 0;
    num_5 = 0;
    %ͳ�ƹ��ж��ٸ�
    for i=1:N
        if sensor_r(1,i) == 7
            num_7 = num_7 + 1;
        elseif sensor_r(1,i) == 6
            num_6 = num_6 + 1;
        else
            num_5 = num_5 + 1;
        end
    end
    
    num_7_sensor = (1:1:num_7);%��һ���� ��Ӧ�����   
    num_6_sensor = ((num_7+1):1:(num_7 + num_6));
    num_5_sensor = ((num_7 + num_6 + 1):1:N);

    matrix_dis_first_any = zeros(N,N);%���ڴ�first��any���ڵ�֮��ľ���
    %��þ��� һ���ڵ� �ֱ���һֻ�ǵ���ÿ���ڵ�ľ���  �����������ͬһ����𣬾�������Ϊinf
    for i=1:N
        for j=1:N
            if ((any(num_7_sensor == i))==1) && ((any(num_7_sensor == j))==1)
             dis_squ = (first_init_wolf(1,i)-any_indivi(1,j))^2 +...
            (first_init_wolf(2,i)-any_indivi(2,j))^2;
            matrix_dis_first_any(i,j) = dis_squ^0.5;
            elseif ((any(num_7_sensor == i))==1) && ((any(num_7_sensor == j))~=1)
               matrix_dis_first_any(i,j) = inf;
            end
            
            if ((any(num_6_sensor == i))==1) && ((any(num_6_sensor == j))==1)
             dis_squ = (first_init_wolf(1,i)-any_indivi(1,j))^2 +...
            (first_init_wolf(2,i)-any_indivi(2,j))^2;
            matrix_dis_first_any(i,j) = dis_squ^0.5;
            elseif ((any(num_6_sensor == i))==1) && ((any(num_6_sensor == j))~=1)
               matrix_dis_first_any(i,j) = inf; 
            end
            
            
            if ((any(num_5_sensor == i))==1) && ((any(num_5_sensor == j))==1)
             dis_squ = (first_init_wolf(1,i)-any_indivi(1,j))^2 +...
            (first_init_wolf(2,i)-any_indivi(2,j))^2;
            matrix_dis_first_any(i,j) = dis_squ^0.5;
            elseif ((any(num_5_sensor == i))==1) && ((any(num_5_sensor == j))~=1)
               matrix_dis_first_any(i,j) = inf; 
            end
            
        end
    end

    % save matrix_dis_first_best.mat matrix_dis_first_best;
    match = zeros(N,2);%��һ�д���ʼ���ڶ��д�Ŀ��λ��
    match(:,1) = (1:1:N);%�����һ�е�˳������
    
    %disp(matrix_dis_first_any);
    [match(:,2),value] = lapjv(matrix_dis_first_any);%����lapjv�㷨
%     disp('���ķѣ�');
%     disp(value);
%     disp('Դͷ    Ŀ��');
%     for i=1:N
%         disp([num2str(match(i,1)),'       ',num2str(match(i,2))]);
%     end
% disp(matrix_dis_first_any);
end
