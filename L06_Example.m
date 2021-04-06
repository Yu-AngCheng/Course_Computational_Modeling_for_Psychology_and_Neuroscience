clear all
close all
dataP = [0.75 ,0.67 ,0.54 ,0.4 ,0.4 ,0.37 ,0.58 ,0.71;
    0.92 ,0.81 ,0.53 ,0.28 ,0.14 ,0.22 ,0.45 ,0.81;
    0.91 ,0.97 ,0.93 ,0.64 ,0.28 ,0.09 ,-0.88 ,0.7;
    0.98 ,0.94 ,0.85 ,0.62 ,0.2 ,0.037 ,0.078 ,0.71;
    0.97 ,0.94 ,0.8 ,0.58 ,0.4 ,0.45 ,0.81 ,0.97;
    0.29 ,0.66 ,0.85 ,0.71 ,0.33 ,0.1 ,0.32 ,0.77];
% number sessions x 10 blocks x 96 trials /(n stimuli)
%训练过程决定了exemplar
Ntrain = ((5*10*96)/8);%训练的时候每个刺激做了600个trial
pfeedback = [ .6 .6 1 1 0 0 .6 .6];%刺激设置的条件
Afeedback = pfeedback.* Ntrain;% 控制的时候用的是频率
feedback= [Afeedback; Ntrain-Afeedback];%第一行是A，第二行是B

Ntest = ((3*10*96)/8);%测试的时候每个刺激做了360个trail
N = repmat(Ntest,1 ,8);
dataF = ceil(Ntest.*(dataP));%测试的时候被试的数据
stimval = linspace (.0625, .9375, 8);%luminance的水平，有8个
%%
% 3 models: GCM DEM GRT; 6 subject
tic
%GCM model
figure;
for i=1:6
    [parameterEst,minusloglikelihood] = ...
        fminsearchbnd(@(c)L06_GCMmodel(c, stimval, dataF(i,:),Ntest, feedback), rand, 0, Inf);
    GCM_c(i)=parameterEst;
    AIC(i,1)=minusloglikelihood+1;
end
for i=1:6
    [~,p]=L06_GCMmodel(GCM_c(i), stimval, dataF(i,:),Ntest, feedback);
    subplot(2,3,i);
    plot(stimval,p,'--+');
    hold on;
    plot(stimval,dataP(i,:),'-*');
    hold off;
    GCM_c_samples=zeros(1,1000);
    parfor j=1:1000
        data_for_bootstrap=binornd(Ntest,p);
        parameterEst = ...
            fminsearchbnd(@(c)L06_GCMmodel(c, stimval, data_for_bootstrap,Ntest, feedback), rand, 0, Inf);
        GCM_c_samples(j)=parameterEst;
    end
    GCM_c_sd(i)=std(GCM_c_samples);
end
% DEM model
figure;
for i=1:6
    parameterEst = ...
        fminsearchbnd(@(param)L06_DEMmodel(param,stimval,dataF(i,:),Ntest, feedback), [rand,rand], [0,0], [Inf,Inf]);
    DEM_c(i)=parameterEst(1);DEM_gamma(i)=parameterEst(2);
    AIC(i,2)=minusloglikelihood+2;
end
for i=1:6
    [~,p]=L06_DEMmodel([DEM_c(i),DEM_gamma(i)], stimval, dataF(i,:),Ntest, feedback);
    subplot(2,3,i);
    plot(stimval,p,'--+');
    hold on;
    plot(stimval,dataP(i,:),'-*');
    hold off;
    DEM_c_samples=zeros(1,1000);DEM_gamma_samples=zeros(1,1000);
    parfor j=1:1000
        data_for_bootstrap=binornd(Ntest,p);
        parameterEst = ...
            fminsearchbnd(@(param)L06_DEMmodel(param,stimval,data_for_bootstrap,Ntest, feedback), [rand,rand], [0,0], [Inf,Inf]);
        DEM_c_samples(j)=parameterEst(1);DEM_gamma_samples(j)=parameterEst(2);
    end
    DEM_c_sd(i)=std(DEM_c_samples);DEM_gamma_sd(i)=std(DEM_gamma_samples);
end
% GRT model
figure;
options=optimset('MaxFunEvals',2000,'MaxIter',2000);
for i=1:6
    parameterEst = ...
        fminsearchbnd(@(param)L06_GRTmodel(param,stimval,dataF(i,:),Ntest), [0.3,0.7,1], [0,0,0], [Inf,Inf,Inf],options);
    GRT_beta1(i)=parameterEst(1);GRT_beta2(i)=parameterEst(2);GRT_sigma(i)=parameterEst(3);
    AIC(i,3)=minusloglikelihood+3;
end
for i=1:6
    [~,p]=L06_GRTmodel([GRT_beta1(i),GRT_beta2(i),GRT_sigma(i)], stimval, dataF(i,:),Ntest);
    subplot(2,3,i);
    plot(stimval,p,'--+');
    hold on;
    plot(stimval,dataP(i,:),'-*');
    hold off;
    parfor j=1:1000
        data_for_bootstrap=binornd(Ntest,p);
        parameterEst = ...
            fminsearchbnd(@(param)L06_GRTmodel(param,stimval,data_for_bootstrap,Ntest), [0.3,0.7,1], [0,0,0], [Inf,Inf,Inf],options);
         GRT_beta1_samples(j)=parameterEst(1);GRT_beta2_samples(j)=parameterEst(2);GRT_sigma_samples(j)=parameterEst(3);
    end
    GRT_beta1_sd(i)=std(GRT_beta1_samples);GRT_beta2_sd(i)=std(GRT_beta2_samples);GRT_sigma_sd(i)=std(GRT_sigma_samples);
end
[alpha,exp_r,xp,pxp,bor] = spm_BMS (AIC, 1e6, 1, 0, 1);