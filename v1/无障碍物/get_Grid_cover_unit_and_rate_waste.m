%%�õ�������
function [cover_rate,waste_rate] = get_Grid_cover_unit_and_rate_waste(sensor_mat1,sersors_r,per_sersons_num,per_sersons_radius_type)
    %�˷ѵı��� = �����븲�� - ��ʵ���� - �ص����ǣ�/ ���븲��
    global M;
    global N;
    global Grid_cen_x_and_y;
    global L;
    global W;

    Grid_cover_unit = zeros(1,M);%1*M�����ڱ���2500����������ϸ���
    Grid_cover_bool = zeros(M,N);%ÿ���������ļ��ĸ���  ֻ����0����
    for i=1:M
        for j=1:N
            if ((Grid_cen_x_and_y(1,i)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(2,i)-sensor_mat1(2,j))^2)...
                    <=sersors_r(1,j)^2
                Grid_cover_bool(i,j) = 1;%�����i�������ĵ���Ա���j���������ڵ㸲��
            end
        end
    end

    %%�������Ϸֲ�����
    for i=1:M
        p = 1;%����Ϊ1
        for j=1:N
            p = p*(1-Grid_cover_bool(i,j));%���ݹ�ʽ����
        end
        Grid_cover_unit(1,i) = 1-p;%���ϸ��Ǹ���
    end

    %���㸲����
    p_sum = 0;%�ۼ���������ĸ��Ǹ���
    for i=1:M
        p_sum = p_sum + Grid_cover_unit(1,i);
    end;
    cover_rate = p_sum*(L*W/M)/(L*W);%������
    
    
    %��ʵ�ĸ������
    area_real = p_sum;
    %disp(['��ʵ�����',num2str(area_real)]);
    
    %ͳ���ص�����
    area_over = 0;%�ص����
    Grid_cover_over_matri = zeros(M,1);
    for i=1:M
        Grid_cover_over_matri(i,1) = sum(Grid_cover_bool(i,:));%ͳ��ÿ�б����Ƕ��ٴ�
    end
    %����ÿ�����㱻�ظ����ǵĴ���
    for i=1:M
        if Grid_cover_over_matri(i,1) >= 2
            area_over  = area_over + Grid_cover_over_matri(i,1) - 1;%�ظ� n-1��
        end
    end
    
    %disp(['�ص������',num2str(area_over)]);
    
    
    %��������������㣬��Ȼ�п��ܳ��ָ��� ��Ϊ�������Щ���
    r_max_area = 158;
    r_mid_area = 117;
    r_min_area = 81;
    %�������ĸ������
    area_ideal = r_max_area*per_sersons_num(1,1) + r_mid_area*per_sersons_num(1,2) + r_min_area*per_sersons_num(1,3);
    
    

    %disp(['���������',num2str(area_ideal)]);
    waste_rate = (area_ideal - area_real - area_over)/area_ideal;
end