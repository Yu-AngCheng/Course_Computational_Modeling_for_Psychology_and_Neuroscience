clear all
close all
dataP = [0.75 ,0.67 ,0.54 ,0.4 ,0.4 ,0.37 ,0.58 ,0.71;
    0.92 ,0.81 ,0.53 ,0.28 ,0.14 ,0.22 ,0.45 ,0.81;
    0.91 ,0.97 ,0.93 ,0.64 ,0.28 ,0.09 ,0.12 ,0.7;
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
%% Maximum likelihood estimation
for modelToFit={'GCM','GRT','DEM'}
    figure;
    for ppt=1:6
        switch char(modelToFit)
            case 'GCM'
                f=@(pars) DEMlnL([pars,1], stimval,...
                    feedback, dataF (ppt,:) , N);
                [theta(ppt,:) ,lnL(ppt) ,exitflag(ppt)]...
                    =fminbnd(f, 0, 100);
            case 'GRT'
                f=@(pars) GRTlnL(pars, stimval,...
                    dataF(ppt,:), N);
                [theta(ppt,:) ,lnL(ppt) ,exitflag(ppt)]...
                    =fminsearchbnd ( f , [0.3,0.7,1] ,...
                    [-1,-1,eps], [2,2,10]);
            case 'DEM'
                f=@(pars) DEMlnL(pars, stimval,...
                    feedback, dataF(ppt,:) , N);
                [theta(ppt,:) ,lnL(ppt) ,exitflag(ppt)]...
                    =fminsearchbnd( f, [5,1] , [0,0], [Inf,Inf]);
            otherwise
                error ( 'Unknown model');
        end
        [~, predP(ppt,:)] = f(theta(ppt,:));
        %         hess = hessian(f,theta(ppt,:),le^-3);
        %         hess = hessian(f,theta(ppt,:));
        %         cov = inv(hess);
        %         thetaSE(ppt,:) = sqrt(diag(cov));
        pptLab = { ' SB ' , 'SEH ' , 'VB ' , 'BG ' , 'NV ' , 'LT ' };
        subplot(2,3,ppt);
        plot (stimval, dataP(ppt,:) , '-+ ' );
        hold all
        plot (stimval, predP(ppt,:) , ' --* ' );
        ylim([0,1]);
        xlabel ('Luminance ' );
        ylabel ( 'P(A) ' );
        title(char(pptLab{ppt}));
        set(gcf, 'Name', char(modelToFit));
%         switch char(modelToFit)
%             case 'GCM'
%                 GCM.theta=theta;
%                 GCM.nlnL=lnL;
%             case 'GRT'
%                 GRT.theta=theta;
%                 GRT.nlnL=lnL;
%             case 'DEM'
%                 DEM.theta=theta;
%                 DEM.nlnL=lnL;
%         end
        clear theta;
    end
end
% for ppt=1:6
%     [[AIC(ppt,:),BIC(ppt,:),AICd(ppt,:)],...
%         [BICd(ppt,:),AICw(ppt,:),BICW(ppt,:)]]=...
%         aicbic([GCM.nlnL(ppt),GRT.nlnL(ppt),...
%         DEM.nlnL(ppt)], [1 3 2],...
%         repmat(Ntest*8,1 ,3));
% end