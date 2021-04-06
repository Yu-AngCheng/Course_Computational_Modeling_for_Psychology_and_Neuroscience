% prime-mask SOA
% metacontrast(mask-target SOA=0) & difference>c
% using E(X)ornot-->important
% based on birth-death-process
clear;
lambdab=0.005;
lambda=1.0;
niu=0.02;
c=48;
blocks=4;
results=[];
EXornot=1;
for cases=1:blocks
    trialeach=25;
    triallist=genTrials(trialeach,[2,2,8]);%1-left,2-right
    trials=length(triallist);
    triallist=[triallist,nan(length(triallist),2)];
    triallist(:,3)=triallist(:,3)*14;
    for i=1:trials
        IL=zeros(1,2);
        ML=zeros(1,2);
        eps=0.1;
        SOAcycles=triallist(i,3)*10;
        ITIcycles=8000;
        IL=zeros(1,2);
        flag=false;
        for j=1:SOAcycles
            IL(triallist(i,1))=1;
            if triallist(i,1)==1
                p1=[ML(1)*niu*eps,(lambda+lambdab)*eps,0];
                p1(3)=1-p1(1)-p1(2);
                if(EXornot)
                    ML(1)=ML(1)+sum(p1.*[-1,1,0]);
                else
                    ML(1)=ML(1)+randsample([-1,1,0],1,true,p1);
                end
                p2=[ML(2)*niu*eps,lambdab*eps,0];
                p2(3)=1-p2(1)-p2(2);
                if(EXornot)
                    ML(2)=ML(2)+sum(p2.*[-1,1,0]);
                else
                    ML(2)=ML(2)+randsample([-1,1,0],1,true,p2);
                end
            elseif triallist(i,1)==2
                p2=[ML(2)*niu*eps,(lambda+lambdab)*eps,0];
                p2(3)=1-p2(1)-p2(2);
                if(EXornot)
                    ML(2)=ML(2)+sum(p2.*[-1,1,0]);
                else
                    ML(2)=ML(2)+randsample([-1,1,0],1,true,p2);
                end
                p1=[ML(1)*niu*eps,lambdab*eps,0];
                p1(3)=1-p1(1)-p1(2);
                if(EXornot)
                    ML(1)=ML(1)+sum(p1.*[-1,1,0]);
                else
                    ML(1)=ML(1)+randsample([-1,1,0],1,true,p1);
                end
            end
        end
        IL=zeros(1,2);
        for j=1:ITIcycles
            IL(triallist(i,2))=1;
            if triallist(i,2)==1
                p1=[ML(1)*niu*eps,(lambda+lambdab)*eps,0];
                p1(3)=1-p1(1)-p1(2);
                if(EXornot)
                    ML(1)=ML(1)+sum(p1.*[-1,1,0]);
                else
                    ML(1)=ML(1)+randsample([-1,1,0],1,true,p1);
                end
                p2=[ML(2)*niu*eps,lambdab*eps,0];
                p2(3)=1-p2(1)-p2(2);
               if(EXornot)
                    ML(2)=ML(2)+sum(p2.*[-1,1,0]);
                else
                    ML(2)=ML(2)+randsample([-1,1,0],1,true,p2);
                end
            elseif triallist(i,2)==2
                p2=[ML(2)*niu*eps,(lambda+lambdab)*eps,0];
                p2(3)=1-p2(1)-p2(2);
                if(EXornot)
                    ML(2)=ML(2)+sum(p2.*[-1,1,0]);
                else
                    ML(2)=ML(2)+randsample([-1,1,0],1,true,p2);
                end
                p1=[ML(1)*niu*eps,lambdab*eps,0];
                p1(3)=1-p1(1)-p1(2);
                if(EXornot)
                    ML(1)=ML(1)+sum(p1.*[-1,1,0]);
                else
                    ML(1)=ML(1)+randsample([-1,1,0],1,true,p1);
                end
            end
            if(ML(1)-ML(2)>c)&&(flag==false)
                triallist(i,4)=j;
                triallist(i,5)=1;
                flag=true;
                break
            elseif(ML(2)-ML(1)>c)&&(flag==false)
                triallist(i,4)=j;
                triallist(i,5)=2;
                flag=true;
                break
            end
        end
    end
    results=[results;triallist];
end
results(:,6)=results(:,1)==results(:,2);
idx=results(:,5)==results(:,2);
[meanresult,sem,group]=grpstats(results(idx,4),results(idx,[6,3]),{'nanmean','sem','gname'});
group=str2double(group);
data=(-meanresult([9,10,11,12,13,14,15,16])+meanresult([1,2,3,4,5,6,7,8]))*0.1;
figure

plot(14:14:112,meanresult(1:8),'r-o');
hold on
plot(14:14:112,meanresult(9:16),'b-o');
box off;

yyaxis right;
plot(14:14:112,data,'k-o');
box off;
legend({'Incongruent','Congruent','Priming effects'},'Location','bestoutside','orientation','horizontal');
legend boxoff;