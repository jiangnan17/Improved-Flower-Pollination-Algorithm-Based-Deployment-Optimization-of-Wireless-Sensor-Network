%һֱһ��ֱ�ߣ����ֱ��ƽ�е�ֱ�ߣ��Ҿ������Ϊr�й�ʽ�ģ����������������������ѡ��һ����ߵ�
% D=|b-c|/��(1+k^2)
% c=b��D��(1+k^2),
% ��һ��ƽ���ߵĺ�����ϵʽy=kx+b��D��(1+k^2)
function [k,b] = get_new_function(k,b,case_b,any_indivi_r)
    if case_b == 1
        b = b - any_indivi_r*rand(1,1) * (1+k^2)^0.5;%�����1�����
    else
        b = b + any_indivi_r *rand(1,1)* (1+k^2)^0.5;%���
    end
end