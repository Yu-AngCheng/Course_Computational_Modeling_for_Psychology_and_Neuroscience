%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
lambda=5;
theta=exp(lambda);
samplesize=20;
resample=10000;
rng default; % only for reproductivity
data=random('Poisson',lambda,samplesize,1);
%% parametric
lambda_estimate=mean(data);
rng shuffle;data_bootstrp=random('Poisson',lambda_estimate,samplesize,resample);
theta_distri_para=exp(mean(data_bootstrp,1));
Bounds=quantile(theta_distri_para,[0.025 0.975]);
LB=Bounds(1);
UB=Bounds(2);

% figure
% histogram(theta_distri_para,40)
% hold on
% h1=xline(LB,'Color',[1 0 0],'LineWidth',2);
% h2=xline(UB,'Color',[1 0 0],'LineWidth',2);
% hold off
% box off
% title('Parametric Bootstrp');
% xlabel('\theta');
% ylabel('Frequency');
% legend([h1 h2],{"x="+num2str(LB),"x="+num2str(UB)});
% legend boxoff;
% print('ͼ1-1','-dpng','-r600');
%% non-parametric
data_bootstrp=data;
rng shuffle;[theta_distri_nonpara,samples] = bootstrp(resample,@(x)exp(nanmean(x)),data_bootstrp);
Bounds=quantile(theta_distri_nonpara,[0.025 0.975]);
LB=Bounds(1);
UB=Bounds(2);

figure
histogram(theta_distri_nonpara,40)
hold on
h1=xline(LB,'Color',[1 0 0],'LineWidth',2);
h2=xline(UB,'Color',[1 0 0],'LineWidth',2);
hold off
box off
title('Non-parametric Bootstrp');
xlabel('\theta');
ylabel('Frequency');
legend([h1 h2],{"x="+num2str(LB),"x="+num2str(UB)});
legend boxoff;
% print('ͼ1-2','-dpng','-r600');
%% true
rng shuffle;data_bootstrp=random('Poisson',lambda,samplesize,resample);
theta_distri=exp(mean(data_bootstrp,1));
Bounds=quantile(theta_distri,[0.025 0.975]);
LB=Bounds(1);
UB=Bounds(2);

figure
histogram(theta_distri,40)
hold on
h1=xline(LB,'Color',[1 0 0],'LineWidth',2);
h2=xline(UB,'Color',[1 0 0],'LineWidth',2);
hold off
box off
title('True Distribution');
xlabel('\theta');
ylabel('Frequency');
legend([h1 h2],{"x="+num2str(LB),"x="+num2str(UB)});
legend boxoff;
% print('ͼ1-3','-dpng','-r600');
