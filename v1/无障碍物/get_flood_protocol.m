function num_send_receive =  get_flood_protocol(any_indivi,sersor_r)
 [~,sersor_size] = size(any_indivi);
 %E = 3*E(send) + 2*E(receive) + 1(固定噪声能耗)
 %构造一个这样的数组 第一列为节点的序号  第二列为发送  第三列为接收  第四列为固定噪声能耗
 num_send_receive = zeros(sersor_size,4);
 num_send_receive(:,1) = (1:1:sersor_size);%初始化列
 num_send_receive(:,2) = 1;%发送个数
 num_send_receive(:,4) = 1;%每个结点固定的消耗能量
 connect_num = zeros(sersor_size,2);%两列，第一列对应节点的序号，第二列对应邻居节点的个数
 connect_num(:,1) = (1:1:sersor_size);%初始化列
 %注意这里有坑，我们求其接收的，相当于求其邻居个数，接收的相当于是其邻居可以发送的。
 for i=1:sersor_size
    for j=1:sersor_size%与除了本身之外的结点进行求距离
        if(i==j)
            continue;
        end
        %这里应该有改进  a->b  b->a  通信距离 是否覆盖是不一样  现统一为 ri + rj
        if ((any_indivi(1,i) - any_indivi(1,j))^2 + (any_indivi(2,i) - any_indivi(2,j))^2)...
                ^0.5 <= 2*sersor_r(1,j);
            connect_num(i,2) = connect_num(i,2) + 1;
        end
    end
 end
    
    num_send_receive(:,3) = connect_num(:,2); %得到邻居节点的个数
end