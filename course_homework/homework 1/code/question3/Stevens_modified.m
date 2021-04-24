function f = Stevens_modified(beta,I, data)
alpha = sum((I.^beta).*data)/sum((I.^beta).^2);
f = sum((data-alpha*I.^beta).^2);
end