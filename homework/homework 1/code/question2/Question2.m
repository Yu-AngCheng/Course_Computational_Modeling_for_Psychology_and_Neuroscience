%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
intensity_level=1:100;
beta_true=0.1;
alpha_true=45;
p_intensity=1./(1+exp(-beta_true.*(intensity_level-alpha_true)));
rng default;% only for reproductivity
data=binornd(1,p_intensity);
%%
figure
plot(intensity_level,p_intensity,'LineWidth',1.5);
hold on;
plot(intensity_level,data,'ro','MarkerFaceColor','r','MarkerEdgeColor','k');
hold off
box off
xlabel('Intensity Level');
ylabel('Probability Correct');
legend({'psychometric function','simulated data'});
legend boxoff
% print('ͼ2-1','-dpng','-r600');
%%
LB = [-Inf, -Inf];
UB = [Inf, Inf];
x0 = [rand, rand];
paramsEst = fminsearchbnd(@(params)mll(params, intensity_level, data), x0, LB, UB);
alphaHat = paramsEst(1);
betaHat = paramsEst(2);
p_intensity_Est=1./(1+exp(-betaHat.*(intensity_level-alphaHat)));

figure
plot(intensity_level,p_intensity_Est,'LineWidth',1.5);
hold on;
plot(intensity_level,p_intensity,'Color','c','LineWidth',1.5);
plot(intensity_level,data,'o','MarkerFaceColor','r','MarkerEdgeColor','k');
hold off
box off
xlabel('Intensity Level');
ylabel('Probability Correct');
legend({'True','Estimated','simulated data'});
legend boxoff
% print('ͼ2-2','-dpng','-r600');
%% parametric bootstrapping
resample=1000;
paramsHat_distri_para=zeros(resample,2);%column 1 for alpha 2 for beta
for i=1:resample
    rng shuffle;data_bootstrp=binornd(1,p_intensity_Est);
    paramsHat_distri_para(i,:)=paraestimation(data_bootstrp);
end
Bounds=quantile(paramsHat_distri_para,[0.025 0.975]);
LB_alphaHat=Bounds(1,1);
UB_alphaHat=Bounds(2,1);
LB_betaHat=Bounds(1,2);
UB_betaHat=Bounds(2,2);