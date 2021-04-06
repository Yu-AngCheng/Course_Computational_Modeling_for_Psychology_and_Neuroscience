function minusloglikelihood = Fechner_psychfun(params,I, data)
alpha = params(1);
beta = params(2);
sigema=params(3);
for i=1:numel(I)
    p(i)=normpdf(data(i),alpha*log(I(i)+beta),sigema);
end
p(p<1e-16) = 1e-16;%为了避免函数在0或1处log出现问题
minusloglikelihood = -sum(log(p));

% minusloglikelihood = 0.5*sum(((data-alpha.*log(I+beta))./sigema).^2)+numel(data)...
% *log(sqrt(2*pi)*sigema);
end
