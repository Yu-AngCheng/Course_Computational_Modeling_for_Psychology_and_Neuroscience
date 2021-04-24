%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 2e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 1; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('nx',nx,'nz',nz,'k',k);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.x = round(rand(1,nx));
    S.z = round(rand(1,nz));
    S.alpha = 1/2;
    S.beta = 1/2;
    init0(i) = S;
end
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, 'TwoCountryQuiz_Question2.txt'), ...
    init0, ...
    'doparallel' , doparallel, ...
    'nchains', nchains,...
    'nburnin', nburnin,...
    'nsamples', nsamples, ...
    'thin', nthin, ...
    'monitorparams', {'x','z','alpha','beta'}, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'cleanup' , 0);
%%
%plot each person's posterior probability
probability=zeros(8,1);
for i=1:8
    data=samples.x(:,:,i);
    probability(i)=sum(data(:))/(nchains*nsamples);
end
figure
xx=1:8;
idx=probability>0.5;
bar(xx(idx),probability(idx),'FaceColor',[1,0.3,0.3],'EdgeColor',[1,0.1,0.1]);
hold on
bar(xx(~idx),probability(~idx),'FaceColor',[0.3,0.3,1],'EdgeColor',[0.1,0.1,1]);
box off;
xticks(xx);
ylim([0 1.2]);
xlabel('No. of people','FontSize',14); 
ylabel('Posterior Probability of being Moldovan','FontSize',14);
%plot each question's posterior probability
probability=zeros(8,1);
for i=1:8
    data=samples.z(:,:,i);
    probability(i)=sum(data(:))/(nchains*nsamples);
end
figure
xx=1:8;
idx=probability>0.5;
bar(xx(idx),probability(idx),'FaceColor',[1,0.3,0.3],'EdgeColor',[1,0.1,0.1]);
hold on
bar(xx(~idx),probability(~idx),'FaceColor',[0.3,0.3,1],'EdgeColor',[0.1,0.1,1]);
box off;
xticks(xx);
ylim([0 1.2]);
xlabel('No. of Question','FontSize',14); 
ylabel('Posterior Probability of being Moldovan','FontSize',14);