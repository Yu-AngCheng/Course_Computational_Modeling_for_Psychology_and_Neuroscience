function minusloglikelihood = mll(theta,n,data)
    len=length(data);
    p=0;
    for i=15:25
        p=p+nchoosek(n,i)*theta.^i*(1-theta).^(n-i);
    end
    p=repmat(p,len-1,1);
    p(len)=nchoosek(n,data(len))*theta.^data(len)*(1-theta).^(n-data(len));
    minusloglikelihood=-sum(log(p));
end

