%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
k = [1 0 0 1 1 0 0 1;
    1 0 0 1 1 0 0 1;
    0 1 1 0 0 1 0 0;
    0 1 1 0 0 1 1 0;
    1 0 0 1 1 0 0 1;
    0 0 0 1 1 0 0 1;
    0 1 0 0 0 1 1 0;
    0 1 1 1 0 1 1 0];
% Constants
[nx,nz] = size(k); % Number of people and questions
nchains = 10; % How Many Chains?
nburnin = 1e5; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 2; % How Often is a Sample Recorded?
doparallel = 1; % Parallel Option
ndic=1;

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('nx',nx,'nz',nz,'k',k);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.x = round(rand(1,nx));
    S.z = round(rand(1,nz));
    S.alpha = 1/2;
    S.beta1 = 1/2;
    S.beta2 = 1/2;
    init0(i) = S;
end
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, 'Copy_of_TwoCountryQuiz_Question3.txt'), ...
    init0, ...
    'doparallel' , doparallel, ...
    'nchains', nchains,...
    'nburnin', nburnin,...
    'nsamples', nsamples, ...
    'thin', nthin, ...
    'monitorparams', {'x','z','alpha','beta1','beta2'}, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'dic', ndic, ... 
    'cleanup' , 0);
