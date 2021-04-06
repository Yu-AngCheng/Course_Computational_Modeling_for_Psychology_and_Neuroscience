function [minusloglikelihood,p] = L06_GCMmodel(c,stimuluslevel,data,N,exemplar)
for i=1:length(stimuluslevel)
    d=abs(stimuluslevel(i)-stimuluslevel);
    s=exp(-c*d);
    sa=sum(s.*exemplar(1,:));
    sb=sum(s.*exemplar(2,:));
    p(i)=sa/(sa+sb);
end
minusloglikelihood=-sum(data.*log(p)+(N-data).*log(1-p));
end

