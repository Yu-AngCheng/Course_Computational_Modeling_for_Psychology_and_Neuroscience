% DEM on Kruschke Data

clear;

sampler = 1; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
load KruschkeData y d1 d2 n nstim nsubj a x;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('y',y,'nstim',nstim,'t',n*nsubj,'a',a,'d1',d1,'d2',d2);

% Initial Values to Supply to WinBugs
for i=1:nchains
	S.w = 0.5;
	S.c = 1;
    S.temp=1;
	init0(i) = S;
end

if ~run_model
    load GCM_1 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
[samples, stats] = matbugs(datastruct, ...
	fullfile(pwd, 'GCM_1.txt'), ...
	'init', init0, ...
	'nChains', nchains, ...
	'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
	'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
	'monitorParams', {'c','w','predy'}, ...
	'Bugdir', 'C:/Program Files/WinBUGS14');
 toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'DEM_1.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams', {'c','w','predy','pr','gamma'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'ndic' ,1);
        toc
    end
%% calculate the DIC and WAIC

    for i=1:n
        temp=samples.pr(:,:,i);
        probability(i,:)=temp(:);
    end
    meandeviance=-2*mean(sum(log(probability)));% also equals to mean(samples.deviance(:));
    pdic=2*var(sum(log(probability)));% also equals to var(samples.deviance(:));
    DIC=meandeviance+pdic;
    lppd=sum(log(mean(probability,2)));
%     pWAIC=sum(var(log(probability),[],2));
    pWAIC=2*sum(log(mean(probability,2))-mean(log(probability),2));
    WAIC=-2*lppd+2*pWAIC;
%     save DEM_1
end

%%
% Posterior Over Gamma
% gamma_distri = sort(reshape(samples.gamma,1,[]));
% Bounds=quantile(gamma_distri,[0.025 0.975]);
% LB=Bounds(1);
% UB=Bounds(2);
% figure;hold on;
% eps=0.1; binsc = 0+eps/2:eps:100-eps/2; binse = 0:eps:100;
% count = histcounts(gamma_distri,binse);
% count = count/sum(count)/eps;
% ph = plot(binsc,count,'k-');
% [gamma,ind] = max(count);
% th = text(25,gamma*1.1,sprintf('%1.3f - %1.3f',LB,UB));
% set(th,'hor','cen','fontsize',14);
% th = text(25,gamma*1.2,sprintf('%d%%',95));
% set(th,'hor','cen','fontsize',14);
% set(gca,'ylim',get(gca,'ylim')*1.4);
% xlabel('\gamma','fontsize',16);
% ylabel('Posterior Density','fontsize',16);
% hold off
% % Predictions
% figure;hold on;
% set(gcf,'units','norm','pos',[.2 .2 .6 .45],'paperpositionmode','auto');
% axis([.5 nstim+.5 -.5 n+.5]);
% xlabel('Stimulus','fontsize',16);
% ylabel('Category Decision','fontsize',16);
% set(gca,'ticklength',[0 0],'fontsize',14,'xtick',1:nstim,'ytick',[0 n],'yticklabel',{'B','A'},'box','off');
% plot(1:nstim,y/nsubj,'k-+');
% errorbar(1:nstim,stats.mean.predy/nsubj,...
% (stats.ci_low.predy-stats.mean.predy)/nsubj,...
% (stats.ci_high.predy-stats.mean.predy)/nsubj,'k:*','LineWidth',1);
% legend({'Raw data','DEM predictions'});
% legend boxoff;
% hold off;