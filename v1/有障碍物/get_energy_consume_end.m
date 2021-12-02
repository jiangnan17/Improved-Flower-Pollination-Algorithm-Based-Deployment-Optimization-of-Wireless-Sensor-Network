%得到消耗的能量
%1.移动的距离的消耗  1mw/m
%2.接收发送的消耗  接收2、发送3mw/连接
%3.辐射面积的消耗   这是固定的 一个圈圈为1mw/m2
%表示优化部署后的能量消耗。
function [match,energy,energy_rate,dis_match_first_best,energy_send_receive_match,energy_precep_match] =  get_energy_consume_end(first_init_wolf,best_indivi,sersor_r,energy_init)
    
    energy_ideal = sum(energy_init);%获得初值能量之和
    [match,value,matrix_dis_first_any] = get_match_value(first_init_wolf,best_indivi,sersor_r);
    %得到指派的所有距离
    distance = value;%得到总的消耗距离
    energy1 = distance * 0.0002;
    %disp(['移动距离的消耗：',num2str(energy1)]);
    [~,N] = size(first_init_wolf);%得到大小
    dis_match_first_best = zeros(N,2);%存匹配后各点之间距离
    dis_match_first_best(:,1) = (1:1:N);%第一列为顺序  用于存匹配后每个节点移动的距离
    for i=1:N
        dis_match_first_best(i,2) = matrix_dis_first_any(match(i,1),match(i,2));
    end
%     disp(dis_match_first_any);
%     pause(100);
   %E = 3*E(send) + 2*E(receive) + 1(固定噪声能耗)
   %构造一个这样的数组 第一列为节点的序号  第二列为发送  第三列为接收  第四列为固定噪声能耗
   num_send_receive =  get_flood_protocol(best_indivi,sersor_r);%
   
   energy_send_receive_match = zeros(N,2);%匹配后的发送、接收的能量消耗
   energy_send_receive_match(:,1) = num_send_receive(:,1);
   for i=1:N
       energy_send_receive_match(i,2) = 0.0005*sum(num_send_receive(i,2)) + 0.0005*sum(num_send_receive(i,3)) + 0.00005* sum(num_send_receive(i,4));
   end
   
   
   %进行还原匹配
   energy_send_receive_match_temp = energy_send_receive_match;
   for i=1:N
       energy_send_receive_match(i,2) = energy_send_receive_match_temp(match(i,2),2);
   end
   
   energy2 = 0.0005*sum(num_send_receive(:,2)) + 0.0005*sum(num_send_receive(:,3)) + 0.00005* sum(num_send_receive(:,4));
   % disp(['接收发送的消耗：',num2str(energy2)]);
   
   [area_radio,radio_area_arr]  = get_precep_energy(sersor_r);
   energy3 =area_radio * 0.0003;%辐射/感知消耗的能量
    %disp(['辐射的消耗：',num2str(energy3)]);
    
   %匹配好后的感知能耗，其实都不用匹配，因为都为1
   energy_precep_match = zeros(N,2);
   energy_precep_match(:,1) = (1:1:N);
   energy_precep_match(:,2) = radio_area_arr;
    


    energy = energy1 + energy2 + energy3;
    %disp(['总的消耗：',num2str(energy)]);
    
    energy_rate = energy/energy_ideal;
    %disp(['能量消耗率：',num2str(energy_rate)]);
end