%一直一个直线，求该直线平行的直线，且距离相差为r有公式的，但求出来的是两条，所以选择一条左边的
% D=|b-c|/√(1+k^2)
% c=b±D√(1+k^2),
% 另一条平行线的函数关系式y=kx+b±D√(1+k^2)
function [k,b] = get_new_function(k,b,case_b,any_indivi_r)
    if case_b == 1
        b = b - any_indivi_r*rand(1,1) * (1+k^2)^0.5;%如果是1则相减
    else
        b = b + any_indivi_r *rand(1,1)* (1+k^2)^0.5;%相加
    end
end