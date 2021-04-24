%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
x=1:10;
samplesize=50;
intensity_level=10;
Stevens.alpha=2;
Stevens.beta=2;
Stevens.sigema=0.4;
rng shuffle
Stevens.error=normrnd(0,Stevens.sigema,samplesize,intensity_level);
Stevens.y=Stevens.alpha*(x.^Stevens.beta)+Stevens.error;
Stevens.aev=mean(Stevens.y);

Fechner.alpha=2;
Fechner.beta=4;
Fechner.sigema=0.4;
rng shuffle
Fechner.error=normrnd(0,Fechner.sigema,samplesize,intensity_level);
Fechner.y=Fechner.alpha.*log(x+Fechner.beta)+Fechner.error;
Fechner.aev=mean(Fechner.y);

figure
plot(x,Stevens.aev,'-o','LineWidth',1.5);
hold on
plot(x,Fechner.aev,'-o','LineWidth',1.5);
hold off
box off
legend({"Stevens' model","Fechner's model"});
legend boxoff
% print('图3-1','-dpng','-r600');
%%
for i=1:50
    data=Stevens.y(i,:);
    LB = 0;
    UB = Inf;
    
    x0 = 2;
    beta = fminsearchbnd(@(beta) Stevens_modified(beta,x,data), x0, LB, UB);
    %Stevens_modified 为仅含beta的一元函数
    %当minus loglikelihood取最小值时，Stevens_modified也取最小值
    %下同
    alpha = sum((x.^beta).*data)/sum((x.^beta).^2);
    sigema=sqrt(sum((data-alpha*(x.^beta)).^2)/numel(data));
    minusloglikelihood=Stevens_psychfun([alpha,beta,sigema],x, data);
    Stevens.fitStevens(i,1)=minusloglikelihood;
    Stevens.fitStevens(i,2)=2*minusloglikelihood+3*log(10);
    
    x0 = 4;
    beta = fminsearchbnd(@(beta) Fechner_modified(beta,x,data), x0, LB, UB);
    %Fechner_modified 为仅含beta的一元函数
    %当minus loglikelihood取最小值时，Fechner_modified也取最小值
    %下同
    alpha = sum(log(x+beta).*data)/sum(log(x+beta).^2);
    sigema=sqrt(sum((data-alpha*log(x+beta)).^2)/numel(data));
    minusloglikelihood=Fechner_psychfun([alpha,beta,sigema],x, data);
    Stevens.fitFechner(i,1)=minusloglikelihood;
    Stevens.fitFechner(i,2)=2*minusloglikelihood+3*log(10);
end
for i=1:50
    data=Fechner.y(i,:);
    LB = 0;
    UB = Inf;
    
    x0 = 2;
    beta = fminsearchbnd(@(beta) Stevens_modified(beta,x,data), x0, LB, UB);
    alpha = sum((x.^beta).*data)/sum((x.^beta).^2);
    sigema=sqrt(sum((data-alpha*(x.^beta)).^2)/numel(data));
    minusloglikelihood=Stevens_psychfun([alpha,beta,sigema],x, data);
    Fechner.fitStevens(i,1)=minusloglikelihood;
    Fechner.fitStevens(i,2)=2*minusloglikelihood+3*log(10);
    
    x0 = 4;
    beta = fminsearchbnd(@(beta) Fechner_modified(beta,x,data), x0, LB, UB);
    alpha = sum(log(x+beta).*data)/sum(log(x+beta).^2);
    sigema=sqrt(sum((data-alpha*log(x+beta)).^2)/numel(data));
    minusloglikelihood=Fechner_psychfun([alpha,beta,sigema],x, data);
    Fechner.fitFechner(i,1)=minusloglikelihood;
    Fechner.fitFechner(i,2)=2*minusloglikelihood+3*log(10);
end

% Data from Stevens
temp1=mean(Stevens.fitStevens);
temp2=mean(Stevens.fitFechner);
temp3=sum(Stevens.fitStevens(:,2)<Stevens.fitFechner(:,2))/50; %Stevens中Stevens'model更好的百分比
% Data from Fechner's
temp4=mean(Fechner.fitStevens);
temp5=mean(Fechner.fitFechner);
temp6=sum(Fechner.fitFechner(:,2)<Fechner.fitStevens(:,2))/50; %Fechner中Fechner's model更好的百分比

%only to organize data
answer1=table([[temp1(1);temp4(1)],[temp2(1);temp5(1)]],[[temp1(2);temp4(2)],[temp2(2);temp5(2)]],...
[[temp3;1-temp6],[1-temp3;temp6]],'VariableNames',{'Mean_minuslnL','Mean_BIC','Percentage_of_best_fit'},...
 'RowNames',{'Data_from_Stevens','Data_from_Fechner'});
%% 
Stevens.CVStevens=zeros(50,1);Stevens.CVFechner=zeros(50,1);
Fechner.CVFechner=zeros(50,1);Fechner.CVStevens=zeros(50,1);
for i=1:50
    data=Stevens.y(i,:);
    CV=0;
    for j=1:10 %data(j)被用作test，其余用来train，下同
        data_CV=[data(1:j-1),data(j+1:10)];
        x=[1:j-1,j+1:10];
        LB = 0;
        UB = Inf;
       
        x0 = 2;
        beta = fminsearchbnd(@(beta) Stevens_modified(beta,x,data_CV), x0, LB, UB);
        alpha = sum((x.^beta).*data_CV)/sum((x.^beta).^2);
        sigema=sqrt(sum((data_CV-alpha*(x.^beta)).^2)/numel(data_CV));
        paramsEst=[alpha,beta,sigema];
        Stevens.CVStevens(i)= Stevens.CVStevens(i)+Stevens_psychfun(paramsEst,j,data(j));
        
        x0 = 4;
        LB = 0;
        UB = Inf;
        beta = fminsearchbnd(@(beta) Fechner_modified(beta,x,data_CV), x0, LB, UB);
        alpha = sum(log(x+beta).*data_CV)/sum(log(x+beta).^2);
        sigema=sqrt(sum((data_CV-alpha*log(x+beta)).^2)/numel(data_CV));
        paramsEst=[alpha,beta,sigema];
        Stevens.CVFechner(i)=Stevens.CVFechner(i)+Fechner_psychfun(paramsEst,j,data(j));
    end
end
for i=1:50
    data=Fechner.y(i,:);
    CV=0;
    for j=1:10
        data_CV=[data(1:j-1),data(j+1:10)];
        x=[1:j-1,j+1:10];
        LB = 0;
        UB = Inf;
        
        x0 = 2;
        beta = fminsearchbnd(@(beta) Stevens_modified(beta,x,data_CV), x0, LB, UB);
        alpha = sum((x.^beta).*data_CV)/sum((x.^beta).^2);
        sigema=sqrt(sum((data_CV-alpha*(x.^beta)).^2)/numel(data_CV));
        paramsEst=[alpha,beta,sigema];
        Fechner.CVStevens(i)= Fechner.CVStevens(i)+Stevens_psychfun(paramsEst,j,data(j));
        
        x0 = 4;
        beta = fminsearchbnd(@(beta) Fechner_modified(beta,x,data_CV), x0, LB, UB);
        alpha = sum(log(x+beta).*data_CV)/sum(log(x+beta).^2);
        sigema=sqrt(sum((data_CV-alpha*log(x+beta)).^2)/numel(data_CV));
        paramsEst=[alpha,beta,sigema];
        Fechner.CVFechner(i)=Fechner.CVFechner(i)+Fechner_psychfun(paramsEst,j,data(j));
    end
end

% Data from Stevens
temp1=mean(Stevens.CVStevens);
temp2=mean(Stevens.CVFechner);
temp3=sum(Stevens.CVStevens<Stevens.CVFechner)/50; %Stevens中Stevens'model更好的百分比
% Data from Fechner's
temp4=mean(Fechner.CVStevens);
temp5=mean(Fechner.CVFechner);
temp6=sum(Fechner.CVFechner<Fechner.CVStevens)/50; %Fechner中Fechner's model更好的百分比

% only to organize data
answer2=table([[temp1;temp4],[temp2;temp5]],[[temp3;1-temp6],[1-temp3;temp6]],...
'VariableNames',{'Mean_CV','Percentage_of_best_fit'},...
 'RowNames',{'Data_from_Stevens','Data_from_Fechner'});