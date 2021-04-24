function minusloglikelihood = mll(params,I,C)
alpha=params(1);
beta=params(2);
pc=1./(1+exp(-beta.*(I-alpha)));
pc(pc<1e-16) = 1e-16;%Ϊ�˱��⺯����0��1��log��������
pc(pc>1-1e-16) = 1-1e-16;%Ϊ�˱��⺯����0��1��log��������
minusloglikelihood=-sum(C.*log(pc)+(1-C).*log(1-pc));
end

