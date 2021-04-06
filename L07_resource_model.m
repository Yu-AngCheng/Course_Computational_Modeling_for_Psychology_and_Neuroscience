clear
traillist=genTrials(1000,360); % 360 angles and each for 1000 repetitions
sd=deg2rad(22);% target error 
K=1/(sd*sd);% transformation to precision
beta=0.15;% weight for binding error
gamma=0.15;% weight for miss error
% pdf of von Mises distribution, a pdf similar to normal distribution but
% particular for those on a ring.
pdf=@(x,miu,precision) exp(precision.*cos(x-miu))./(2*pi*besseli(0,precision));

eps=0.01;% steps
theta=0:eps:2*pi;

for i=1:length(traillist)
    target=deg2rad(traillist(i,1));
    % 3 distractors and 1 target
    distractor1=rand()*2*pi;
    distractor2=rand()*2*pi;
    distractor3=rand()*2*pi;
    % The first item is for target error, the second item it for missing error
    % and the third is for binding error
    p=(1-gamma-beta)*pdf(theta,target,K)+ ... 
        gamma*1/(2*pi)+...
        beta/3*(pdf(theta,distractor1,K)+pdf(theta,distractor2,K)+pdf(theta,distractor3,K));
    traillist(i,2)=randsample(rad2deg(theta),1,true,p);
    % de-centered
    traillist(i,3)=traillist(i,2)-distractor1/pi*180;
    traillist(i,4)=traillist(i,2)-distractor2/pi*180;
    traillist(i,5)=traillist(i,2)-distractor3/pi*180;
end

% Distance from target
DistancefromT=traillist(:,2)-traillist(:,1);
% wrap up for those larget than 180 or less than -180
DistancefromT(DistancefromT>180)=DistancefromT(DistancefromT>180)-360;
DistancefromT(DistancefromT<-180)=DistancefromT(DistancefromT<-180)+360;
% draw a picture
figure
histogram(DistancefromT(:),100);

% Distance from Distractors
DistancefromNT=traillist(:,[3,4,5]);
% wrap up for those larget than 180 or less than -180
idx_1=DistancefromNT>180;
DistancefromNT(idx_1)=DistancefromNT(idx_1)-360;
idx_2=DistancefromNT<-180;
DistancefromNT(idx_2)=DistancefromNT(idx_2)+360;
% draw a picture
figure
histogram(DistancefromNT(:),100);
xlim([-180,180]);