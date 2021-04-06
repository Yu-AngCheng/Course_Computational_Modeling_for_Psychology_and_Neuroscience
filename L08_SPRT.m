% Sequential Probability Ratio Test
miu1=4;
miu2=2;
sigma1=3;
sigma2=3;
% H0 x from population 1
% H1 x from population 2
log_prior_odds=log(1);
log_posterior_odds=log_prior_odds;
for i=1:200000
    if(log_posterior_odds>30)
        fprintf("H1 is correct\n");
        break;
    elseif(log_posterior_odds<-30)
        fprintf("H0 is correct\n");
        break;
    end
    x=normrnd(miu2,sigma2);
    log_likelihood_ratio=log(normpdf(x,miu2,sigma2)/normpdf(x,miu1,sigma1));
    log_posterior_odds=log_posterior_odds+log_likelihood_ratio;
end