function Y=logistic(X,AL)
K=4.52;
g=AL(3)*K;
Y=1./(1+exp(-X.*g));
end