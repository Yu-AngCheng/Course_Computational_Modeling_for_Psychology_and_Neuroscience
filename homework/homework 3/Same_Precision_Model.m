%% The Seven Scientists

clear;

%% Data
x = [-27.020,3.570,8.191,9.898,9.603,9.945,10.056];

% Constants
n = length(x);

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e5;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 1; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    S.mu = 0; % An Intial Value
    S.lambda =1; % Intial Values For All The Precisions
    init0(i) = S;
end

% Use JAGS to Sample
tic
fprintf( 'Running JAGS ...\n' );
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, 'SevenScientists_Question2.txt'), ...
    init0, ...
    'doparallel' , doparallel, ...
    'nchains', nchains,...
    'nburnin', nburnin,...
    'nsamples', nsamples, ...
    'thin', nthin, ...
    'monitorparams', {'mu','sigma','pr'}, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'cleanup' , 0 , ...
    'dic' ,1);
toc
% calculate the DIC and WAIC
DIC=stats.dic;
for i=1:n
    temp=samples.pr(:,:,i);
    probability(i,:)=temp(:);
end
lppd=sum(log(mean(probability,2)));
pWAIC=sum(var(log(probability),[],2));
WAIC=-2*lppd+2*pWAIC;
%%
% Posterior Over Rate
miu_distri = sort(reshape(samples.mu,1,[]));
Bounds=quantile(miu_distri,[0.025 0.975]);
LB=Bounds(1);
UB=Bounds(2);
figure;hold on;
eps=0.5; binsc = -30+eps/2:eps:40-eps/2; binse = -30:eps:40;
count = histcounts(miu_distri,binse);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
[miu,ind] = max(count);
th = text(0.2,miu*1.1,sprintf('%1.3f - %1.3f',LB,UB));
set(th,'hor','cen','fontsize',14);
th = text(0.2,miu*1.2,sprintf('%d%%',95));
set(th,'hor','cen','fontsize',14);
set(gca,'ylim',get(gca,'ylim')*1.4);
xlabel('\mu','fontsize',16);
ylabel('Posterior Density','fontsize',16);
title('Same Measurement Precision Model','fontsize',16);
hold off
save same_measurement_precision_model