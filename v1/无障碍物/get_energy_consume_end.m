%�õ����ĵ�����
%1.�ƶ��ľ��������  1mw/m
%2.���շ��͵�����  ����2������3mw/����
%3.�������������   ���ǹ̶��� һ��ȦȦΪ1mw/m2
%��ʾ�Ż��������������ġ�
function [match,energy,energy_rate,dis_match_first_best,energy_send_receive_match,energy_precep_match] =  get_energy_consume_end(first_init_wolf,best_indivi,sersor_r,energy_init)
    
    energy_ideal = sum(energy_init);%��ó�ֵ����֮��
    [match,value,matrix_dis_first_any] = get_match_value(first_init_wolf,best_indivi,sersor_r);
    %�õ�ָ�ɵ����о���
    distance = value;%�õ��ܵ����ľ���
    energy1 = distance * 0.0002;
    %disp(['�ƶ���������ģ�',num2str(energy1)]);
    [~,N] = size(first_init_wolf);%�õ���С
    dis_match_first_best = zeros(N,2);%��ƥ������֮�����
    dis_match_first_best(:,1) = (1:1:N);%��һ��Ϊ˳��  ���ڴ�ƥ���ÿ���ڵ��ƶ��ľ���
    for i=1:N
        dis_match_first_best(i,2) = matrix_dis_first_any(match(i,1),match(i,2));
    end
%     disp(dis_match_first_any);
%     pause(100);
   %E = 3*E(send) + 2*E(receive) + 1(�̶������ܺ�)
   %����һ������������ ��һ��Ϊ�ڵ�����  �ڶ���Ϊ����  ������Ϊ����  ������Ϊ�̶������ܺ�
   num_send_receive =  get_flood_protocol(best_indivi,sersor_r);%
   
   energy_send_receive_match = zeros(N,2);%ƥ���ķ��͡����յ���������
   energy_send_receive_match(:,1) = num_send_receive(:,1);
   for i=1:N
       energy_send_receive_match(i,2) = 0.0005*sum(num_send_receive(i,2)) + 0.0005*sum(num_send_receive(i,3)) + 0.00005* sum(num_send_receive(i,4));
   end
   
   
   %���л�ԭƥ��
   energy_send_receive_match_temp = energy_send_receive_match;
   for i=1:N
       energy_send_receive_match(i,2) = energy_send_receive_match_temp(match(i,2),2);
   end
   
   energy2 = 0.0005*sum(num_send_receive(:,2)) + 0.0005*sum(num_send_receive(:,3)) + 0.00005* sum(num_send_receive(:,4));
   % disp(['���շ��͵����ģ�',num2str(energy2)]);
   
   [area_radio,radio_area_arr]  = get_precep_energy(sersor_r);
   energy3 =area_radio * 0.0003;%����/��֪���ĵ�����
    %disp(['��������ģ�',num2str(energy3)]);
    
   %ƥ��ú�ĸ�֪�ܺģ���ʵ������ƥ�䣬��Ϊ��Ϊ1
   energy_precep_match = zeros(N,2);
   energy_precep_match(:,1) = (1:1:N);
   energy_precep_match(:,2) = radio_area_arr;
    


    energy = energy1 + energy2 + energy3;
    %disp(['�ܵ����ģ�',num2str(energy)]);
    
    energy_rate = energy/energy_ideal;
    %disp(['���������ʣ�',num2str(energy_rate)]);
end