%spike_train_simulation
clear
% stimulus
eps=0.001;
t=eps:eps:0.5;
stimulus=[repmat(-40,1,length(t)/5),repmat(-20,1,length(t)/5),...
    repmat(0,1,length(t)/5),repmat(20,1,length(t)/5),...
    repmat(40,1,length(t)/5)]; %#ok<*REPMAT>
% tuning function
rmax=55;
smax=0;
sigmaf=20;
tuning_function=@(s) rmax*exp(-0.5*(s-smax).^2/sigmaf^2);
% mean firing rate
firing_rate=tuning_function(stimulus);
% spike_train
p=firing_rate*eps;
num=5000;
p=repmat(p,num,1);
spike_train=binornd(1,p,num,length(t));
% figure
figure;
plot(t,mean(spike_train)/eps);
hold on;
plot(t,firing_rate,'LineWidth',2);
hold off;