function [match,value,matrix_dis_first_any] = get_match_value(first_init_wolf,any_indivi,sensor_r)

    [~,N] = size(first_init_wolf);%得到有多少个节点
    
    num_7 = 0;%半径为7的个数
    num_6 = 0;
    num_5 = 0;
    %统计共有多少个
    for i=1:N
        if sensor_r(1,i) == 7
            num_7 = num_7 + 1;
        elseif sensor_r(1,i) == 6
            num_6 = num_6 + 1;
        else
            num_5 = num_5 + 1;
        end
    end
    
    num_7_sensor = (1:1:num_7);%第一种类 对应的序号   
    num_6_sensor = ((num_7+1):1:(num_7 + num_6));
    num_5_sensor = ((num_7 + num_6 + 1):1:N);

    matrix_dis_first_any = zeros(N,N);%用于存first到any各节点之间的距离
    %算好距离 一个节点 分别到另一只狼当中每个节点的距离  如果不是属于同一个类别，距离设置为inf
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
    match = zeros(N,2);%第一列存起始，第二列存目的位置
    match(:,1) = (1:1:N);%解决第一列的顺序问题
    
    %disp(matrix_dis_first_any);
    [match(:,2),value] = lapjv(matrix_dis_first_any);%调用lapjv算法
%     disp('共耗费：');
%     disp(value);
%     disp('源头    目的');
%     for i=1:N
%         disp([num2str(match(i,1)),'       ',num2str(match(i,2))]);
%     end
% disp(matrix_dis_first_any);
end
