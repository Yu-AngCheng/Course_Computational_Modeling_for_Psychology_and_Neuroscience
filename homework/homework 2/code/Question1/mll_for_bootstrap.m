function minusloglikelihood = mll_for_bootstrap(theta,n,data)
    p=data.*log(theta)+(n-data).*log(1-theta);
    minusloglikelihood=-sum(p);
end
