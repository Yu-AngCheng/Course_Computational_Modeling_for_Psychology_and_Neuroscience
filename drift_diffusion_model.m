%% DDM
clear
a = 1; % upper bound
z = 0.5;% percentage of starting point
e = 0;% end
drift_rate=0.5;
steps=0.001;
s = 1;% standard error of one step
current = a*z;% real starting point
time = 0;

t(1) = time;
evidence(1) = current;
n = 2;% just for index
while (current<a)&&(current>e)
    time=time+steps;
    current=current+normrnd(drift_rate*steps,sqrt(s*s*steps));
    t(n)=time;
    evidence(n)=current;
    n=n+1;
end
figure
plot(t,evidence);
ylim([e,a]);
