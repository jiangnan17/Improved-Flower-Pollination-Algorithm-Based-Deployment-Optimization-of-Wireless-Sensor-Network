%得到消耗的能量
%1.移动的距离的消耗  1mw/m
%2.接收发送的消耗  接收2、发送3mw/连接
%3.辐射面积的消耗   这是固定的 一个圈圈为1mw/m2
function [energy,energy_rate] =  get_energy_consume(first_init_wolf,any_indivi,serson_r,energy_init)
    
    energy_ideal = sum(energy_init);%获得初值能量之和
    [~,value] = get_match_value(first_init_wolf,any_indivi,serson_r);
    %得到指派的所有距离
    distance = value;%得到总的消耗距离
    energy1 = distance * 0.0002;
    %disp(['移动距离的消耗：',num2str(energy1)]);
    
   %E = 3*E(send) + 2*E(receive) + 1(固定噪声能耗)
   %构造一个这样的数组 第一列为节点的序号  第二列为发送  第三列为接收  第四列为固定噪声能耗
   num_send_receive =  get_flood_protocol(any_indivi,serson_r);%best_indivi是一个整的结构体类型
   energy2 = 0.0005*sum(num_send_receive(:,2)) + 0.0005*sum(num_send_receive(:,3)) + 0.00005* sum(num_send_receive(:,4));
   % disp(['接收发送的消耗：',num2str(energy2)]);
    
   [area_radio,~]  = get_precep_energy(serson_r);
   energy3 = area_radio * 0.0003;%辐射/感知消耗的能量
    %disp(['辐射的消耗：',num2str(energy3)]);
    
    


    energy = energy1 + energy2 + energy3;
    %disp(['总的消耗：',num2str(energy)]);
    
    energy_rate = energy/energy_ideal;
    %disp(['能量消耗率：',num2str(energy_rate)]);
end