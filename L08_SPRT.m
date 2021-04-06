%% Sequential Probability Ratio Test
% parameters for 2 normal distributions
miu_1=4;
miu_2=2;
sigma_1=3;
sigma_2=3;

% H0: sample x is from population 1
% H1: sample x is from population 2
log_prior_odds=log(1);
log_posterior_odds=log_prior_odds;
for i=1:200000
    % set the criteria
    if(log_posterior_odds>30)
        fprintf("H1 is correct\n");
        break;
    elseif(log_posterior_odds<-30)
        fprintf("H0 is correct\n");
        break;
    end
    % get a sample
    x=normrnd(miu_2,sigma_2);
    % calculate the likelihood
    log_likelihood_ratio=log(normpdf(x,miu_2,sigma_2)/normpdf(x,miu_1,sigma_1));
    % update the posterior
    log_posterior_odds=log_posterior_odds+log_likelihood_ratio;
end