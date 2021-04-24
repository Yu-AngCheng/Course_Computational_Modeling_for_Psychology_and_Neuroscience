function f = Fechner_modified(beta,I, data)
alpha = sum(log(I+beta).*data)/sum(log(I+beta).^2);
f = sum((data-log(I+beta)*alpha).^2);
end