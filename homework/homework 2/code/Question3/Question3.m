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
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
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
    fullfile(pwd, 'TwoCountryQuiz_Question3.txt'), ...
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

beta1_distri = sort(reshape(samples.beta1,1,[]));
Bounds=quantile(beta1_distri,[0.025 0.975]);
LB_theta=Bounds(1);
UB_theta=Bounds(2);
figure;hold on;
eps=.004; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
count = histc(beta1_distri,binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
[beta1,ind] = max(count);
th = text(0.2,beta1*1.1,sprintf('%1.3f - %1.3f',LB_theta,UB_theta));
set(th,'hor','cen','fontsize',14);
th = text(0.2,beta1*1.2,sprintf('%d%%',95));
set(th,'hor','cen','fontsize',14);
set(gca,'ylim',get(gca,'ylim')*1.4);
xlabel('\beta_{1}','fontsize',16);
ylabel('Posterior Density','fontsize',16);
hold off
beta2_distri = sort(reshape(samples.beta2,1,[]));
Bounds=quantile(beta2_distri,[0.025 0.975]);
LB_theta=Bounds(1);
UB_theta=Bounds(2);

figure;hold on;
eps=.004; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
count = histc(beta2_distri,binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
[beta2,ind] = max(count);
th = text(0.2,beta2*1.1,sprintf('%1.3f - %1.3f',LB_theta,UB_theta));
set(th,'hor','cen','fontsize',14);
th = text(0.2,beta2*1.2,sprintf('%d%%',95));
set(th,'hor','cen','fontsize',14);
set(gca,'ylim',get(gca,'ylim')*1.4);
xlabel('\beta_{2}','fontsize',16);
ylabel('Posterior Density','fontsize',16);
hold off
