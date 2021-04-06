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
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    S.mu = 0; % An Intial Value
    S.phig = 0.5*ones(1,n);
    S.lambdatemp = 1; % Intial Values For All The Precisions
    S.delta=1;
%     S.zg = floor(rand(1,n)*2);
    init0(i) = S;
end
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'SevenScientists_low_high_precision.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'mu','phig','zg','sigma','pr'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'ndic',1,...
        'workingdir' , 'tmpjags');
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
save Low_High_Measurement_Precision_Model1
