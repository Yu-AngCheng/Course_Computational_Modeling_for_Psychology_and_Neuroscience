%% DDM
clear
upper = 1; % upper bound
start_percent = 0.5;% percentage of starting point
lower = 0;% end
drift_rate=0.5;
steps=0.001;
s = 1;% standard error of one step
current = upper*start_percent;% real starting point
time = 0;

t(1) = time;
evidence(1) = current;
while (current<upper)&&(current>lower)
    time=time+steps;
    current=current+normrnd(drift_rate*steps,sqrt(s*s*steps));
    t = [t,time];
    evidence = [evidence,current];
end
figure
plot(t,evidence);
ylim([lower,upper]);
%%
clear
start_percent = 0.5;% percentage of starting point
lower = 0;% end
steps=0.001;
s = 1;% standard error of one step
allruns = 10000;
uppers = 1:0.2:2.8;
drift_rate = 0.1;

for i = 1:10
    upper = uppers(i);
    ACC(i) = 0;
    for run = 1:allruns
        time = 0;
        current = upper*start_percent;% starting point
        while (current<upper)&&(current>lower)
            time=time+steps;
            current=current+normrnd(drift_rate*steps,sqrt(s*s*steps));
        end
        if current>=upper
            ACC(i) = ACC(i) + 1;
        end
    end
end
plot(uppers,ACC/allruns);