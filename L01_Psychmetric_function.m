% 心理物理学曲线的参数拟合
I = [5, 7, 8, 10, 12]';
N = [10, 12, 10, 11, 9]';
C = [2, 3, 7, 9, 9]';
LB = [0, 0];
UB = [Inf, Inf];
x0 = [rand, rand];
[paramsEst, minuslli, exitflag] = ...
fminsearchbnd(@(params)PsychFun(params, I, N, C), x0, LB, UB);
alphaHat = paramsEst(1);
betaHat = paramsEst(2);
lli = - minuslli;


function minuslli = PsychFun(params, I, N, C)
alpha = params(1);
beta = params(2);
pc = 1-exp(-(I/alpha).^beta);
pc(pc<1e-16) = 1e-16;%为了避免函数在0或1处log出现问题
pc(pc>1-1e-16) = 1-1e-16;%为了避免函数在0或1处log出现问题
minuslli = -sum(C.*log(pc)+(N-C).*log(1-pc));
end