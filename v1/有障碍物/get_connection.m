%%�ж��Ƿ���ͨ
function [is_connec,adjacencyMatrix,adjacencyMatrix_dis ]= get_connection(any_indivi,sersors_r)

[~,serson_size] = size(any_indivi);%����õ��������ڵ�Ĵ�С  serson_size����ǰ��N


%%������ͨ��
adjacencyMatrix=zeros(serson_size,serson_size);%���崫���������ڽӾ���  ��������ͼ
adjacencyMatrix_dis=ones(serson_size,serson_size).*inf;%���崫���������ڽӾ���ľ��� ��Ҳ������ͼ,����תΪ����ͼ����
for i=1:1:serson_size
    for j=(i+1):1:serson_size   %��Ϊ�ж��ְ뾶  ��ͨ�ý��д���
        %û���Լ�ָ���Լ���
%         if i == j
%             continue;
%         end
        
%         dis_squ = (any_indivi(1,i)-any_indivi(1,j))^2 +...
%         (any_indivi(2,i)-any_indivi(2,j))^2;
%                     %�Ѿ������
%             if dis_squ <= (2 * sersors_r(1,i))^2;%������֮����Ը�֪
%                 adjacencyMatrix(i,j)=1;%����ͼ   
%                 %�������
%                 adjacencyMatrix_dis(i,j)=dis_squ^0.5;%����ͼ
%             end
            

            flag1 = 0;%��ʾA-->B����ͨ
            flag2 = 0;%��ʾB-->A����ͨ
            dis_squ = (any_indivi(1,i)-any_indivi(1,j))^2 +...
        (any_indivi(2,i)-any_indivi(2,j))^2;
           %�Ѿ������
            if dis_squ <= (2 * sersors_r(1,i))^2;%������֮����Ը�֪
                flag1 = 1;%��A-->B���򵽴�
            end
            if dis_squ <= (2 * sersors_r(1,j))^2;%������֮����Ը�֪
                flag2 = 1;%��B-->A���򵽴�
            end
            
            %���������ͼ����A-->B  B-->A����������ͼ  �պ÷����������󣬵�����������
            if flag1==1 && flag2==1;%������֮����Ի����֪
                adjacencyMatrix(i,j)=1;%����ͼ
                adjacencyMatrix(j,i)=1;
                
                %�������
                adjacencyMatrix_dis(i,j)=dis_squ^0.5;%����ͼ
                adjacencyMatrix_dis(j,i)=dis_squ^0.5;
            end   
    end
end 

%�������һ�ε�����ͼͼ�ڽӾ���
% save adjacencyMatrix_dis.mat adjacencyMatrix_dis;
% save adjacencyMatrix.mat adjacencyMatrix;

S=zeros(serson_size,serson_size);
for m=1:1:serson_size-1
    S=S+adjacencyMatrix^m;  %ȫ���ӵ�S
end;

%S=M+M^2+M^3+...+M^(N-1)��
%����N��M������������
%��S����Ԫ��Ϊ�㣬����ͨ��
%S�����㣬����ͨ��
%�������ڽӾ���ͼ����ͨ���ж�׼�򡷲鵽��
is_connec = 0;%��ʼ��Ϊ����ͨ
if all(all(S))==1 %�ж�S�����ǲ���ȫΪ����Ԫ��
    is_connec = 1;%ȫΪ��0  ����ͨ
end
end