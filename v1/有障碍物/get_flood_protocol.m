function num_send_receive =  get_flood_protocol(any_indivi,sersor_r)
 [~,sersor_size] = size(any_indivi);
 %E = 3*E(send) + 2*E(receive) + 1(�̶������ܺ�)
 %����һ������������ ��һ��Ϊ�ڵ�����  �ڶ���Ϊ����  ������Ϊ����  ������Ϊ�̶������ܺ�
 num_send_receive = zeros(sersor_size,4);
 num_send_receive(:,1) = (1:1:sersor_size);%��ʼ����
 num_send_receive(:,2) = 1;%���͸���
 num_send_receive(:,4) = 1;%ÿ�����̶�����������
 connect_num = zeros(sersor_size,2);%���У���һ�ж�Ӧ�ڵ����ţ��ڶ��ж�Ӧ�ھӽڵ�ĸ���
 connect_num(:,1) = (1:1:sersor_size);%��ʼ����
 %ע�������пӣ�����������յģ��൱�������ھӸ��������յ��൱�������ھӿ��Է��͵ġ�
 for i=1:sersor_size
    for j=1:sersor_size%����˱���֮��Ľ����������
        if(i==j)
            continue;
        end
        %����Ӧ���иĽ�  a->b  b->a  ͨ�ž��� �Ƿ񸲸��ǲ�һ��  ��ͳһΪ ri + rj
        if ((any_indivi(1,i) - any_indivi(1,j))^2 + (any_indivi(2,i) - any_indivi(2,j))^2)...
                ^0.5 <= 2*sersor_r(1,j);
            connect_num(i,2) = connect_num(i,2) + 1;
        end
    end
 end
    
    num_send_receive(:,3) = connect_num(:,2); %�õ��ھӽڵ�ĸ���
end