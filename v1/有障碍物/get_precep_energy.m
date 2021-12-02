function [area_radio,radio_area_arr] =  get_precep_energy(sersor_r)
    [~,sersor_size] = size(sersor_r);%得到一个个体对应的节点数目
    radio_area_arr = zeros(sersor_size,1);%把每一个个体的感知能耗存储进去
    for i=1:sersor_size
        radio_area_arr(i,1) = pi* sersor_r(1,i)^2;
    end
    area_radio = sum(radio_area_arr);%总的辐射面积
end