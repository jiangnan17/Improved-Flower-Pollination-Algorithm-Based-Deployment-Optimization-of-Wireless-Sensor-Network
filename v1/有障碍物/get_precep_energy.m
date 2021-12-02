function [area_radio,radio_area_arr] =  get_precep_energy(sersor_r)
    [~,sersor_size] = size(sersor_r);%�õ�һ�������Ӧ�Ľڵ���Ŀ
    radio_area_arr = zeros(sersor_size,1);%��ÿһ������ĸ�֪�ܺĴ洢��ȥ
    for i=1:sersor_size
        radio_area_arr(i,1) = pi* sersor_r(1,i)^2;
    end
    area_radio = sum(radio_area_arr);%�ܵķ������
end