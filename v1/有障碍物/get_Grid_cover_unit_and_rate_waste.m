%%�õ�������
function [cover_rate,waste_rate] = get_Grid_cover_unit_and_rate_waste(sensor_mat1,sersors_r,per_sersons_num,per_sersons_radius_type)
    %�˷ѵı��� = �����븲�� - ��ʵ���� - �ص����ǣ�/ ���븲��
    global M;
    global N;
    global Grid_cen_x_and_y;
    global L;
    global W;
    
    %��İ취 �ֱ�д����
    %�������ݵĴ���  ͨ����ȡ
    fiftynum = (1:1:50);
    caseone_row1 = [fiftynum(1:1:15),fiftynum(36:1:50)];
    caseone_row2 = [fiftynum(1:1:14),fiftynum(37:1:50)];
    caseone_row3 = [fiftynum(1:1:13),fiftynum(38:1:50)];
    caseone_row4 = [fiftynum(1:1:12),fiftynum(39:1:50)];
    caseone_row5 = [fiftynum(1:1:11),fiftynum(40:1:50)];
    caseone_row6 = [fiftynum(1:1:10),fiftynum(41:1:50)];
    caseone_row7 = [fiftynum(1:1:9),fiftynum(42:1:50)];
    caseone_row8 = [fiftynum(1:1:8),fiftynum(43:1:50)];
    caseone_row9 = [fiftynum(1:1:7),fiftynum(44:1:50)];
    caseone_row10 = [fiftynum(1:1:6),fiftynum(45:1:50)];
    caseone_row11 = [fiftynum(1:1:5),fiftynum(46:1:50)];
    caseone_row12 = [fiftynum(1:1:4),fiftynum(47:1:50)];
    caseone_row13 = [fiftynum(1:1:3),fiftynum(48:1:50)];
    caseone_row14 = [fiftynum(1:1:2),fiftynum(49:1:50)];
    caseone_row15 = [fiftynum(1:1:1),fiftynum(50:1:50)];



    caseone_row36 = [fiftynum(1:1:1),fiftynum(50:1:50)];
    caseone_row37 = [fiftynum(1:1:2),fiftynum(49:1:50)];
    caseone_row38 = [fiftynum(1:1:3),fiftynum(48:1:50)];
    caseone_row39 = [fiftynum(1:1:4),fiftynum(47:1:50)];
    caseone_row40 = [fiftynum(1:1:5),fiftynum(46:1:50)];
    caseone_row41 = [fiftynum(1:1:6),fiftynum(45:1:50)];
    caseone_row42 = [fiftynum(1:1:7),fiftynum(44:1:50)];
    caseone_row43 = [fiftynum(1:1:8),fiftynum(43:1:50)];
    caseone_row44 = [fiftynum(1:1:9),fiftynum(42:1:50)];
    caseone_row45 = [fiftynum(1:1:10),fiftynum(41:1:50)];
    caseone_row46 = [fiftynum(1:1:11),fiftynum(40:1:50)];
    caseone_row47 = [fiftynum(1:1:12),fiftynum(39:1:50)];
    caseone_row48 = [fiftynum(1:1:13),fiftynum(38:1:50)];
    caseone_row49 = [fiftynum(1:1:14),fiftynum(37:1:50)];
    caseone_row50 = [fiftynum(1:1:15),fiftynum(36:1:50)];
    %���εĴ���
    caseone_row16 = [fiftynum(25:1:25),fiftynum(26:1:26)];
    caseone_row17 = [fiftynum(24:1:25),fiftynum(26:1:27)];
    caseone_row18 = [fiftynum(23:1:25),fiftynum(26:1:28)];
    caseone_row19 = [fiftynum(22:1:25),fiftynum(26:1:29)];
    caseone_row20 = [fiftynum(21:1:25),fiftynum(26:1:30)];
    caseone_row21 = [fiftynum(20:1:25),fiftynum(26:1:31)];
    caseone_row22 = [fiftynum(19:1:25),fiftynum(26:1:32)];
    caseone_row23 = [fiftynum(18:1:25),fiftynum(26:1:33)];
    caseone_row24 = [fiftynum(17:1:25),fiftynum(26:1:34)];
    caseone_row25 = [fiftynum(16:1:25),fiftynum(26:1:35)];
    
    
    caseone_row26 = [fiftynum(16:1:25),fiftynum(26:1:35)];
    caseone_row27 = [fiftynum(17:1:25),fiftynum(26:1:34)];
    caseone_row28 = [fiftynum(18:1:25),fiftynum(26:1:33)];
    caseone_row29 = [fiftynum(19:1:25),fiftynum(26:1:32)];
    caseone_row30 = [fiftynum(20:1:25),fiftynum(26:1:31)];
    caseone_row31 = [fiftynum(21:1:25),fiftynum(26:1:30)];
    caseone_row32 = [fiftynum(22:1:25),fiftynum(26:1:29)];
    caseone_row33 = [fiftynum(23:1:25),fiftynum(26:1:28)];
    caseone_row34 = [fiftynum(24:1:25),fiftynum(26:1:27)];
    caseone_row35 = [fiftynum(25:1:25),fiftynum(26:1:26)];
    
    Grid_cover_unit = zeros(L,W);%1*M�����ڱ���2500����������ϸ���
    Grid_cover_bool = zeros(L,W,N);%ÿ���������ļ��ĸ���  ֻ����0����
    %����Ҳ�ÿ���  �漰�����ص�����
    for i=1:L
		for k=1:W
            
            %���ķ��������ĸ�������
            if i == 1 && (any(caseone_row1==k))==1
                continue;
            end
            if i == 2 && (any(caseone_row2==k))==1
                continue;
            end

            if i == 3 && (any(caseone_row3==k))==1
                continue;
            end
            if i == 4 && (any(caseone_row4==k))==1
                continue;
            end
            if i == 5 && (any(caseone_row5==k))==1
                continue;
            end
            if i == 6 && (any(caseone_row6==k))==1
                continue;
            end
            if i == 7 && (any(caseone_row7==k))==1
                continue;
            end
            if i == 8 && (any(caseone_row8==k))==1
                continue;
            end
            if i == 9 && (any(caseone_row9==k))==1
                continue;
            end
            if i == 10 && (any(caseone_row10==k))==1
                continue;
            end
            if i == 11 && (any(caseone_row11==k))==1
                continue;
            end
            if i == 12 && (any(caseone_row12==k))==1
                continue;
            end
            if i == 13 && (any(caseone_row13==k))==1
                continue;
            end
            if i == 14 && (any(caseone_row14==k))==1
                continue;
            end
            if i == 15 && (any(caseone_row15==k))==1
                continue;
            end


            if i == 36 && (any(caseone_row36==k))==1
                continue;
            end
            if i == 37 && (any(caseone_row37==k))==1
                continue;
            end
            if i == 38 && (any(caseone_row38==k))==1
                continue;
            end
            if i == 39 && (any(caseone_row39==k))==1
                continue;
            end
            if i == 40 && (any(caseone_row40==k))==1
                continue;
            end
            if i == 41 && (any(caseone_row41==k))==1
                continue;
            end
            if i == 42 && (any(caseone_row42==k))==1
                continue;
            end
            if i == 43 && (any(caseone_row43==k))==1
                continue;
            end
            if i == 44 && (any(caseone_row44==k))==1
                continue;
            end
            if i == 45 && (any(caseone_row45==k))==1
                continue;
            end
            if i == 46 && (any(caseone_row46==k))==1
                continue;
            end

            if i == 47 && (any(caseone_row47==k))==1
                continue;
            end
            if i == 48 && (any(caseone_row48==k))==1
                continue;
            end
            if i == 49 && (any(caseone_row49==k))==1
                continue;
            end

            if i == 50 && (any(caseone_row50==k))==1
                continue;
            end

            %�����ϰ�����ж�
            if i == 16 && (any(caseone_row16==k))==1
                continue;
            end
            if i == 17 && (any(caseone_row17==k))==1
                continue;
            end

            if i == 18 && (any(caseone_row18==k))==1
                continue;
            end
            if i == 19 && (any(caseone_row19==k))==1
                continue;
            end
            if i == 20 && (any(caseone_row20==k))==1
                continue;
            end
            if i == 21 && (any(caseone_row21==k))==1
                continue;
            end
            if i == 22 && (any(caseone_row22==k))==1
                continue;
            end
            if i == 23 && (any(caseone_row23==k))==1
                continue;
            end
            if i == 24 && (any(caseone_row24==k))==1
                continue;
            end
            if i == 25 && (any(caseone_row25==k))==1
                continue;
            end
            if i == 26 && (any(caseone_row26==k))==1
                continue;
            end
            if i == 27 && (any(caseone_row27==k))==1
                continue;
            end
            if i == 28 && (any(caseone_row28==k))==1
                continue;
            end
            if i == 29 && (any(caseone_row29==k))==1
                continue;
            end
            if i == 30 && (any(caseone_row30==k))==1
                continue;
            end


            if i == 31 && (any(caseone_row31==k))==1
                continue;
            end
            if i == 32 && (any(caseone_row32==k))==1
                continue;
            end
            if i == 33 && (any(caseone_row33==k))==1
                continue;
            end
            if i == 34 && (any(caseone_row34==k))==1
                continue;
            end
            if i == 35 && (any(caseone_row35==k))==1
                continue;
            end
            
            
            
			for j=1:N
				if ((Grid_cen_x_and_y(i,k,1)-sensor_mat1(1,j))^2 + (Grid_cen_x_and_y(i,k,2)-sensor_mat1(2,j))^2)...
						<=sersors_r(1,j)^2
					Grid_cover_bool(i,k,j) = 1;%�����i�������ĵ���Ա���j���������ڵ㸲��
				end
			end
		end
    end
    
    
    %%�������Ϸֲ�����
    for i=1:L
		for k=1:W
            
            %���ķ��������ĸ�������
            if i == 1 && (any(caseone_row1==k))==1
                continue;
            end
            if i == 2 && (any(caseone_row2==k))==1
                continue;
            end

            if i == 3 && (any(caseone_row3==k))==1
                continue;
            end
            if i == 4 && (any(caseone_row4==k))==1
                continue;
            end
            if i == 5 && (any(caseone_row5==k))==1
                continue;
            end
            if i == 6 && (any(caseone_row6==k))==1
                continue;
            end
            if i == 7 && (any(caseone_row7==k))==1
                continue;
            end
            if i == 8 && (any(caseone_row8==k))==1
                continue;
            end
            if i == 9 && (any(caseone_row9==k))==1
                continue;
            end
            if i == 10 && (any(caseone_row10==k))==1
                continue;
            end
            if i == 11 && (any(caseone_row11==k))==1
                continue;
            end
            if i == 12 && (any(caseone_row12==k))==1
                continue;
            end
            if i == 13 && (any(caseone_row13==k))==1
                continue;
            end
            if i == 14 && (any(caseone_row14==k))==1
                continue;
            end
            if i == 15 && (any(caseone_row15==k))==1
                continue;
            end


            if i == 36 && (any(caseone_row36==k))==1
                continue;
            end
            if i == 37 && (any(caseone_row37==k))==1
                continue;
            end
            if i == 38 && (any(caseone_row38==k))==1
                continue;
            end
            if i == 39 && (any(caseone_row39==k))==1
                continue;
            end
            if i == 40 && (any(caseone_row40==k))==1
                continue;
            end
            if i == 41 && (any(caseone_row41==k))==1
                continue;
            end
            if i == 42 && (any(caseone_row42==k))==1
                continue;
            end
            if i == 43 && (any(caseone_row43==k))==1
                continue;
            end
            if i == 44 && (any(caseone_row44==k))==1
                continue;
            end
            if i == 45 && (any(caseone_row45==k))==1
                continue;
            end
            if i == 46 && (any(caseone_row46==k))==1
                continue;
            end

            if i == 47 && (any(caseone_row47==k))==1
                continue;
            end
            if i == 48 && (any(caseone_row48==k))==1
                continue;
            end
            if i == 49 && (any(caseone_row49==k))==1
                continue;
            end

            if i == 50 && (any(caseone_row50==k))==1
                continue;
            end

            %�����ϰ�����ж�
            if i == 16 && (any(caseone_row16==k))==1
                continue;
            end
            if i == 17 && (any(caseone_row17==k))==1
                continue;
            end

            if i == 18 && (any(caseone_row18==k))==1
                continue;
            end
            if i == 19 && (any(caseone_row19==k))==1
                continue;
            end
            if i == 20 && (any(caseone_row20==k))==1
                continue;
            end
            if i == 21 && (any(caseone_row21==k))==1
                continue;
            end
            if i == 22 && (any(caseone_row22==k))==1
                continue;
            end
            if i == 23 && (any(caseone_row23==k))==1
                continue;
            end
            if i == 24 && (any(caseone_row24==k))==1
                continue;
            end
            if i == 25 && (any(caseone_row25==k))==1
                continue;
            end
            if i == 26 && (any(caseone_row26==k))==1
                continue;
            end
            if i == 27 && (any(caseone_row27==k))==1
                continue;
            end
            if i == 28 && (any(caseone_row28==k))==1
                continue;
            end
            if i == 29 && (any(caseone_row29==k))==1
                continue;
            end
            if i == 30 && (any(caseone_row30==k))==1
                continue;
            end


            if i == 31 && (any(caseone_row31==k))==1
                continue;
            end
            if i == 32 && (any(caseone_row32==k))==1
                continue;
            end
            if i == 33 && (any(caseone_row33==k))==1
                continue;
            end
            if i == 34 && (any(caseone_row34==k))==1
                continue;
            end
            if i == 35 && (any(caseone_row35==k))==1
                continue;
            end
            
			p = 1;%����Ϊ1
			for j=1:N
				p = p*(1-Grid_cover_bool(i,k,j));%���ݹ�ʽ����
			end
			Grid_cover_unit(i,k) = 1-p;%���ϸ��Ǹ���
		end
    end

   %���㸲����
    p_sum = 0;%�ۼ���������ĸ��Ǹ���
    for i=1:L
        for j=1:W
            p_sum = p_sum + Grid_cover_unit(i,j);
        end
    end
    
    
    cover_rate = p_sum*(L*W/M)/((L*W)- 480 - 210);%������%���������ʵ����������480    ������450 ���ݵ�������Ϊ512 ���ε㣺220 ʵ�����200 ����210
    
    
    %��ʵ�ĸ������
    area_real = p_sum;
    %disp(['��ʵ�����',num2str(area_real)]);
    
    %ͳ���ص�����
    area_over = 0;%�ص����
    Grid_cover_over_matri = zeros(L,W);
    for i=1:L
        for j=1:W
            Grid_cover_over_matri(i,j) = sum(Grid_cover_bool(i,j,:));%ͳ�Ƹ����㱻���Ƕ��ٴ�
        end
    end
    %����ÿ�����㱻�ظ����ǵĴ���
    for i=1:L
        for j=1:W
            if Grid_cover_over_matri(i,j) >= 2
                area_over  = area_over + Grid_cover_over_matri(i,j) - 1;%�ظ� n-1��
            end
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