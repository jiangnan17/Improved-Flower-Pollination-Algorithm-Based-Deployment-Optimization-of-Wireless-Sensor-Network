%%判断是否连通
function [is_connec,adjacencyMatrix,adjacencyMatrix_dis ]= get_connection(any_indivi,sersors_r)

[~,serson_size] = size(any_indivi);%计算得到传感器节点的大小  serson_size是以前的N


%%计算连通性
adjacencyMatrix=zeros(serson_size,serson_size);%定义传感器互联邻接矩阵  这是有向图
adjacencyMatrix_dis=ones(serson_size,serson_size).*inf;%定义传感器互联邻接矩阵的距离 这也是有向图,试着转为无向图矩阵
for i=1:1:serson_size
    for j=(i+1):1:serson_size   %因为有多种半径  连通得进行处理
        %没有自己指向自己的
%         if i == j
%             continue;
%         end
        
%         dis_squ = (any_indivi(1,i)-any_indivi(1,j))^2 +...
%         (any_indivi(2,i)-any_indivi(2,j))^2;
%                     %把距离加入
%             if dis_squ <= (2 * sersors_r(1,i))^2;%；两点之间可以感知
%                 adjacencyMatrix(i,j)=1;%有向图   
%                 %加入距离
%                 adjacencyMatrix_dis(i,j)=dis_squ^0.5;%有向图
%             end
            

            flag1 = 0;%表示A-->B不连通
            flag2 = 0;%表示B-->A不连通
            dis_squ = (any_indivi(1,i)-any_indivi(1,j))^2 +...
        (any_indivi(2,i)-any_indivi(2,j))^2;
           %把距离加入
            if dis_squ <= (2 * sersors_r(1,i))^2;%；两点之间可以感知
                flag1 = 1;%可A-->B单向到达
            end
            if dis_squ <= (2 * sersors_r(1,j))^2;%；两点之间可以感知
                flag2 = 1;%可B-->A单向到达
            end
            
            %构造出无向图，当A-->B  B-->A当做是无向图  刚好符合我们需求，但不是最完美
            if flag1==1 && flag2==1;%；两点之间可以互相感知
                adjacencyMatrix(i,j)=1;%无向图
                adjacencyMatrix(j,i)=1;
                
                %加入距离
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%无向图
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end   
    end
end 

%保存最后一次的有向图图邻接矩阵
% save adjacencyMatrix_dis.mat adjacencyMatrix_dis;
% save adjacencyMatrix.mat adjacencyMatrix;

S=zeros(serson_size,serson_size);
for m=1:1:serson_size-1
    S=S+adjacencyMatrix^m;  %全部加到S
end;

%S=M+M^2+M^3+...+M^(N-1)；
%其中N是M的行数或列数
%若S中有元素为零，则不连通；
%S中无零，则连通，
%《基于邻接矩阵图的连通性判定准则》查到的
is_connec = 0;%初始化为不连通
if all(all(S))==1 %判断S向量是不是全为非零元素
    is_connec = 1;%全为非0  则连通
end
end