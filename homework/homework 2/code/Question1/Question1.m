%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
nattempts = 950;
nfails = 949;
n=50;
data=[nan(1,nfails),30];
LB = 0;
UB = 1;
x0 = rand;
thetaHat = fminsearchbnd(@(theta)mll(theta, n, data), x0, LB, UB);
resample=10000;
theta_distri=zeros(1,resample);
parfor i=1:resample
    rng shuffle;
    data_bootstrp=binornd(n,thetaHat,1,nattempts);
    theta_distri(i)=fminsearchbnd(@(theta)mll_for_bootstrap...
    (theta, n, data_bootstrp), x0, LB, UB);
end
Bounds=quantile(theta_distri,[0.025 0.975]);
LB_thetaHat=Bounds(1);
UB_thetaHat=Bounds(2);

figure;hold on;
eps=.004; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
count = histc(theta_distri,binse);count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
theta = max(count);
th = text(thetaHat,theta*1.1,sprintf('%1.3f - %1.3f',LB_thetaHat,UB_thetaHat));
set(th,'hor','cen','fontsize',14);
th = text(thetaHat,theta*1.2,sprintf('%d%%',95));
set(th,'hor','cen','fontsize',14);
set(gca,'ylim',get(gca,'ylim')*1.4);
xlabel('Rate','fontsize',16);
ylabel('Distribution','fontsize',16);
hold off;
%%
nattempts = 5;
nfails = 4;
n = 50; % questions
y = [ones(nfails,1);0]; % Indicate which scores are censored
z = [nan*ones(nfails,1);30]'; % All scores except the last are unknown
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 2e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 1; % Parallel Option
ndic=1;
% Assign Matlab Variables to the Observed Nodes
datastruct = struct('nattempts',nattempts,'n',n,'z',z,'y',y);
% Initialize Unobserved Variables
for i=1:nchains
    S.theta = .8;
    init0(i) = S;
end
% Use JAGS to Sample
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, 'grandpa.txt'), ...
    init0, ...
    'doparallel' , doparallel, ...
    'nchains', nchains,...
    'nburnin', nburnin,...
    'nsamples', nsamples, ...
    'thin', nthin, ...
    'dic', ndic, ... 
    'monitorparams', {'theta','z'}, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'cleanup' , 0);
theta_distri = sort(reshape(samples.theta,1,[]));
Bounds=quantile(theta_distri,[0.025 0.975]);
LB_theta=Bounds(1);
UB_theta=Bounds(2);

% Posterior Over Rate
figure;hold on;
eps=.004; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
count = histc(theta_distri,binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
[theta,ind] = max(count);
th = text(binsc(ind),theta*1.1,sprintf('%1.3f - %1.3f',LB_theta,UB_theta));
set(th,'hor','cen','fontsize',14);
th = text(binsc(ind),theta*1.2,sprintf('%d%%',95));
set(th,'hor','cen','fontsize',14);
set(gca,'ylim',get(gca,'ylim')*1.4);
xlabel('Rate','fontsize',16);
ylabel('Posterior Density','fontsize',16);
hold off