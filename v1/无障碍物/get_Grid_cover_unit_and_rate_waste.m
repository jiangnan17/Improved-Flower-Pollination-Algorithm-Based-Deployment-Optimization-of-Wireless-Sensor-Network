%%得到覆盖率
function [cover_rate,waste_rate] = get_Grid_cover_unit_and_rate_waste(sensor_mat1,sersors_r,per_sersons_num,per_sersons_radius_type)
    %浪费的比例 = （理想覆盖 - 真实覆盖 - 重叠覆盖）/ 理想覆盖
    global M;
    global N;
    global Grid_cen_x_and_y;
    global L;
    global W;

    Grid_cover_unit = zeros(1,M);%1*M个用于保存2500个网格的联合概率
    Grid_cover_bool = zeros(M,N);%每个网格中心监测的概率  只能是0或者
    for i=1:M
        for j=1:N
            if ((Grid_cen_x_and_y(1,i)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(2,i)-sensor_mat1(2,j))^2)...
                    <=sersors_r(1,j)^2
                Grid_cover_bool(i,j) = 1;%代表第i网格中心点可以被第j个传感器节点覆盖
            end
        end
    end

    %%计算联合分布概率
    for i=1:M
        p = 1;%设置为1
        for j=1:N
            p = p*(1-Grid_cover_bool(i,j));%根据公式计算
        end
        Grid_cover_unit(1,i) = 1-p;%联合覆盖概率
    end

    %计算覆盖率
    p_sum = 0;%累加所有网格的覆盖概率
    for i=1:M
        p_sum = p_sum + Grid_cover_unit(1,i);
    end;
    cover_rate = p_sum*(L*W/M)/(L*W);%覆盖率
    
    
    %真实的覆盖面积
    area_real = p_sum;
    %disp(['真实面积：',num2str(area_real)]);
    
    %统计重叠区域
    area_over = 0;%重叠面积
    Grid_cover_over_matri = zeros(M,1);
    for i=1:M
        Grid_cover_over_matri(i,1) = sum(Grid_cover_bool(i,:));%统计每行被覆盖多少次
    end
    %计算每个监测点被重复覆盖的次数
    for i=1:M
        if Grid_cover_over_matri(i,1) >= 2
            area_over  = area_over + Grid_cover_over_matri(i,1) - 1;%重复 n-1次
        end
    end
    
    %disp(['重叠面积：',num2str(area_over)]);
    
    
    %利用网格点来计算，不然有可能出现负数 因为网格点有些误差
    r_max_area = 158;
    r_mid_area = 117;
    r_min_area = 81;
    %存放理想的覆盖面积
    area_ideal = r_max_area*per_sersons_num(1,1) + r_mid_area*per_sersons_num(1,2) + r_min_area*per_sersons_num(1,3);
    
    

    %disp(['理想面积：',num2str(area_ideal)]);
    waste_rate = (area_ideal - area_real - area_over)/area_ideal;
end