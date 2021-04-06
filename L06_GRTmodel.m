function [minusloglikelihood,p] = L06_GRTmodel(param,stimuluslevel,data,N)
beta1=param(1);beta2=param(2);sigma=param(3);
if(beta1>=beta2)
    minusloglikelihood=realmax;
    return;
end
p=normcdf((beta1-stimuluslevel)./sigma)+1-normcdf((beta2-stimuluslevel)./sigma);
p(p<1e-16) = 1e-16;
p(p>1-1e-16) = 1-1e-16;
minusloglikelihood=-sum(data.*log(p)+(N-data).*log(1-p));
end

