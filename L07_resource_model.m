clear
traillist=genTrials(1000,360); 
sd=22/360*2*pi;
K=1/(sd*sd);% transformation to von mises distribution
beta=0.15;% weight for binding error
gamma=0.15;% weight for miss
p=1-beta-gamma;% weight for target
eps=0.01;% steps
theta=0:eps:2*pi;
myfun=@(x,miu,precision) exp(precision.*cos(x-miu))./(2*pi*besseli(0,precision));
for i=1:length(traillist)
    target=traillist(i,1)/180*pi;distractor1=rand()*2*pi;
    distractor2=rand()*2*pi;distractor3=rand()*2*pi;
    p=(1-gamma-beta)*myfun(theta,target,K)+gamma*1/(2*pi)+...
    beta/3*(myfun(theta,distractor1,K)+myfun(theta,distractor2,K)+myfun(theta,distractor3,K));
    traillist(i,2)=randsample(theta/(2*pi)*360,1,true,p);
    traillist(i,3)=traillist(i,2)-distractor1/pi*180;traillist(i,4)=traillist(i,2)-distractor2/pi*180;
    traillist(i,5)=traillist(i,2)-distractor3/pi*180;
end
% error=traillist(:,2)-traillist(:,1);
% error(error>180)=error(error>180)-360;
% error(error<-180)=error(error<-180)+360;
% [count,centers]=hist(error,20);
% plot(centers,count,'r-o');
DistancefromNT=traillist(:,[3,4,5]);
idx=DistancefromNT>180;
DistancefromNT(idx)=DistancefromNT(idx)-360;
idxx=DistancefromNT<-180;
DistancefromNT(idxx)=DistancefromNT(idxx)+360;
histogram(DistancefromNT(:),100);
xlim([-180,180]);