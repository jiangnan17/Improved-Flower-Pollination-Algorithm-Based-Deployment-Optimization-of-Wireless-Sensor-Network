function L1 = Levy1(d)
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;%��׼��̬�ֲ���u
    v=randn(1,d);%��׼��̬�ֲ���v
    step=u./abs(v).^(1/beta);
    L1 = 0.001*step;
end