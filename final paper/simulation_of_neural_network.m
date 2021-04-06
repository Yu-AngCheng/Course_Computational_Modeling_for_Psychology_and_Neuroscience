clear
tic
blocks=4;
results=[];
parfor cases=1:blocks
    trialeach=5;
    triallist=genTrials(trialeach,[2,2,160]);%1-left,2-right
    trials=length(triallist);
    triallist=[triallist,nan(trials,2)];
    for i=1:trials
        IL=zeros(1,5);
        RL=zeros(1,5);
        AL=zeros(1,3);
        CL=zeros(1,2);
        ML=zeros(1,2);
        primecycles=43;
%         maskcycles=50;
        maskcycles=100;
        SOAcycles=triallist(i,3)*5;
        targetcycles=200;
        ITIcycles=1000;
        flag=false;
        for j=1:primecycles
            IL(triallist(i,1))=1;
            RL=myupdate(IL,RL,AL,CL,ML,"RL");
            CL=myupdate(IL,RL,AL,CL,ML,"CL");
            ML=myupdate(IL,RL,AL,CL,ML,"ML");
            AL=myupdate(IL,RL,AL,CL,ML,"AL");
        end
%         IL=zeros(1,5);
%         for j=1:triallist(i,3)
%             RL=myupdate(IL,RL,AL,CL,ML,"RL");
%             CL=myupdate(IL,RL,AL,CL,ML,"CL");
%             ML=myupdate(IL,RL,AL,CL,ML,"ML");
%             AL=myupdate(IL,RL,AL,CL,ML,"AL");
%         end
        IL=zeros(1,5);
        for j=1:maskcycles
            IL(3)=1;
            RL=myupdate(IL,RL,AL,CL,ML,"RL");
            CL=myupdate(IL,RL,AL,CL,ML,"CL");
            ML=myupdate(IL,RL,AL,CL,ML,"ML");
            AL=myupdate(IL,RL,AL,CL,ML,"AL");
        end
        IL=zeros(1,5);
        for j=1:SOAcycles
            RL=myupdate(IL,RL,AL,CL,ML,"RL");
            CL=myupdate(IL,RL,AL,CL,ML,"CL");
            ML=myupdate(IL,RL,AL,CL,ML,"ML");
            AL=myupdate(IL,RL,AL,CL,ML,"AL");
        end
        for j=1:targetcycles
            IL(triallist(i,2)+3)=1;
            RL=myupdate(IL,RL,AL,CL,ML,"RL");
            CL=myupdate(IL,RL,AL,CL,ML,"CL");
            ML=myupdate(IL,RL,AL,CL,ML,"ML");
            AL=myupdate(IL,RL,AL,CL,ML,"AL");
            if (ML(1)>0.62)&&(flag==false)
                triallist(i,4)=j;
                flag=true;
                triallist(i,5)=1;
            elseif (ML(2)>0.62)&&(flag==false)
                triallist(i,4)=j;
                flag=true;
                triallist(i,5)=2;
            end
            
        end
        IL=zeros(1,5);
        for j=1:ITIcycles
            RL=myupdate(IL,RL,AL,CL,ML,"RL");
            CL=myupdate(IL,RL,AL,CL,ML,"CL");
            ML=myupdate(IL,RL,AL,CL,ML,"ML");
            AL=myupdate(IL,RL,AL,CL,ML,"AL");
            if (ML(1)>0.62)&&(flag==false)
                triallist(i,4)=j+targetcycles;
                flag=true;
                triallist(i,5)=1;
            elseif (ML(2)>0.62)&&(flag==false)
                triallist(i,4)=j+targetcycles;
                flag=true;
                triallist(i,5)=2;
            end
        end
    end
    results=[results;triallist];
end
toc

results(:,6)=results(:,1)==results(:,2);
idx=results(:,5)==results(:,2);
[meanresult,sem,group]=grpstats(results(:,4),results(:,[6,3]),{'nanmean','sem','gname'});
group=str2double(group);
data=(meanresult(1:160)-meanresult(161:320))*0.1;
figure

plot(5:5:800,meanresult(1:160),'r-o');
hold on
plot(5:5:800,meanresult(161:320),'b-o');
box off;
legend({'Incongruent','Congruent'},'Location','bestoutside','orientation','horizontal');
legend boxoff;
figure
plot(5:5:800,data,'k-o');
box off;
legend({'Priming effects'},'Location','bestoutside','orientation','horizontal');
legend boxoff;