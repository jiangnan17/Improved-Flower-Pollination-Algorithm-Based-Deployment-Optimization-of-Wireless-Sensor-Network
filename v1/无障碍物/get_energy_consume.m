%�õ����ĵ�����
%1.�ƶ��ľ��������  1mw/m
%2.���շ��͵�����  ����2������3mw/����
%3.�������������   ���ǹ̶��� һ��ȦȦΪ1mw/m2
function [energy,energy_rate] =  get_energy_consume(first_init_wolf,any_indivi,serson_r,energy_init)
    
    energy_ideal = sum(energy_init);%��ó�ֵ����֮��
    [~,value] = get_match_value(first_init_wolf,any_indivi,serson_r);
    %�õ�ָ�ɵ����о���
    distance = value;%�õ��ܵ����ľ���
    energy1 = distance * 0.0002;
    %disp(['�ƶ���������ģ�',num2str(energy1)]);
    
   %E = 3*E(send) + 2*E(receive) + 1(�̶������ܺ�)
   %����һ������������ ��һ��Ϊ�ڵ�����  �ڶ���Ϊ����  ������Ϊ����  ������Ϊ�̶������ܺ�
   num_send_receive =  get_flood_protocol(any_indivi,serson_r);%best_indivi��һ�����Ľṹ������
   energy2 = 0.0005*sum(num_send_receive(:,2)) + 0.0005*sum(num_send_receive(:,3)) + 0.00005* sum(num_send_receive(:,4));
   % disp(['���շ��͵����ģ�',num2str(energy2)]);
    
   [area_radio,~]  = get_precep_energy(serson_r);
   energy3 = area_radio * 0.0003;%����/��֪���ĵ�����
    %disp(['��������ģ�',num2str(energy3)]);
    
    


    energy = energy1 + energy2 + energy3;
    %disp(['�ܵ����ģ�',num2str(energy)]);
    
    energy_rate = energy/energy_ideal;
    %disp(['���������ʣ�',num2str(energy_rate)]);
end