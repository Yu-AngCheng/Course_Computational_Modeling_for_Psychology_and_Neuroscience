function [minusloglikelihood,p] = L06_DEMmodel(param,stimuluslevel,data,N,exemplar)
c=param(1);gamma=param(2);
for i=1:length(stimuluslevel)
    d=abs(stimuluslevel(i)-stimuluslevel);
    s=exp(-c*d);
    sa=sum(s.*exemplar(1,:));
    sb=sum(s.*exemplar(2,:));
    p(i)=sa^gamma/(sa^gamma+sb^gamma);
end
p(p<1e-16) = 1e-16;
p(p>1-1e-16) = 1-1e-16;
minusloglikelihood=-sum(data.*log(p)+(N-data).*log(1-p));
end

