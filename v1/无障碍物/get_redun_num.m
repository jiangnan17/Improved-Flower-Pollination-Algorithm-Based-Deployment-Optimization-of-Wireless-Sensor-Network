%%������ڵ����
function get_redun_num()
global N;
global M;
global L;
global W;
global r;
global x;
global y;
global Grid_cover_bool;
global sensor_mat;
global adjacencyMatrix;
global Grid_cover_unit;
global Grid_cen_x_and_y;
global cover_p;
x1 = x;
y1 = y;
%�ֱ�õ�����
Grid_cover_bool1 = Grid_cover_bool;
sensor_mat1 = sensor_mat;
Grid_cen_x_and_y1 = Grid_cen_x_and_y;
Grid_cover_unit1 = Grid_cover_unit;
adjacencyMatrix1 = adjacencyMatrix;

best_sum = 0;
for i=1:1
    sum1 = 0;
    for j=1:M
        Grid_cover_bool1(j,i) = 0;%�Ѹô������ڵ�ȥ�� ����bool�ĸ�֪��Ϊ0
        Grid_cover_unit1 = get_Grid_cover_unit(Grid_cover_unit1);
    end
    p_sum = 0;%�ۼ���������ĸ��Ǹ���
    for i1=1:M
        p_sum = p_sum + Grid_cover_unit(1,i1);
    end;
    cover_p1 = p_sum*10*10/2500;%������
    disp('������Ϊ��');
    disp(cover_p1);
    if cover_p1>cover_p
        sum1 = sum1+1;
    else
        break;
    end
    for z=1:1:N%�����ڽӾ���
        if z~=i
            adjacencyMatrix1(i,z)=0;%����ͼ
            adjacencyMatrix1(z,i)=0;
        end
    end
end
    
    