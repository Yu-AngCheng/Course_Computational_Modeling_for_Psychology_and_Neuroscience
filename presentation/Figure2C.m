clear
%% absolute
for V3temp=0:20:200
    V1=[110,120,130,140];
    V2=[150,150,150,150];
    V3=repmat(V3temp,size(V1));
    V=[V1;V2;V3];
    [firing_rate_mean,noise]=fun(V,"absolute");
    for i=1:4
        p(i)=mvncdf([0,0],...
        [firing_rate_mean(1,i)-firing_rate_mean(2,i),firing_rate_mean(3,i)-firing_rate_mean(2,i)],...
        [noise(1,i)^2+noise(2,i)^2,noise(2,i)^2;noise(2,i)^2,noise(3,i)^2+noise(2,i)^2]);
        q(i)=1-mvncdf([0,0],...
        [firing_rate_mean(1,i)-firing_rate_mean(3,i),firing_rate_mean(2,i)-firing_rate_mean(3,i)],...
        [noise(3,i)^2+noise(1,i)^2,noise(3,i)^2;noise(3,i)^2,noise(3,i)^2+noise(2,i)^2]);
    end
    
    V1=[160,170,180,190];
    V2=[150,150,150,150];
    V3=repmat(V3temp,size(V1));
    V=[V1;V2;V3];
    [firing_rate_mean,noise]=fun(V,"absolute");
    for i=1:4
        p(i+4)=mvncdf([0,0],...
        [firing_rate_mean(2,i)-firing_rate_mean(1,i),firing_rate_mean(3,i)-firing_rate_mean(1,i)],...
        [noise(1,i)^2+noise(2,i)^2,noise(1,i)^2;noise(1,i)^2,noise(3,i)^2+noise(1,i)^2]);
        q(i+4)=1-mvncdf([0,0],...
        [firing_rate_mean(1,i)-firing_rate_mean(3,i),firing_rate_mean(2,i)-firing_rate_mean(3,i)],...
        [noise(3,i)^2+noise(1,i)^2,noise(3,i)^2;noise(3,i)^2,noise(3,i)^2+noise(2,i)^2]);
    end
    answer_absolute(V3/20+1)=mean(p./q);
    clear p q;
end
%% normalization
for V3temp=0:20:200
    V1=[110,120,130,140];
    V2=[150,150,150,150];
    V3=repmat(V3temp,size(V1));
    V=[V1;V2;V3];
    [firing_rate_mean,noise]=fun(V,"relative");
    for i=1:4
        p(i)=mvncdf([0,0],...
        [firing_rate_mean(1,i)-firing_rate_mean(2,i),firing_rate_mean(3,i)-firing_rate_mean(2,i)],...
        [noise(1,i)^2+noise(2,i)^2,noise(2,i)^2;noise(2,i)^2,noise(3,i)^2+noise(2,i)^2]);
        q(i)=1-mvncdf([0,0],...
        [firing_rate_mean(1,i)-firing_rate_mean(3,i),firing_rate_mean(2,i)-firing_rate_mean(3,i)],...
        [noise(3,i)^2+noise(1,i)^2,noise(3,i)^2;noise(3,i)^2,noise(3,i)^2+noise(2,i)^2]);
    end
    
    V1=[160,170,180,190];
    V2=[150,150,150,150];
    V3=repmat(V3temp,size(V1));
    V=[V1;V2;V3];
    [firing_rate_mean,noise]=fun(V,"relative");
    for i=1:4
        p(i+4)=mvncdf([0,0],...
        [firing_rate_mean(2,i)-firing_rate_mean(1,i),firing_rate_mean(3,i)-firing_rate_mean(1,i)],...
        [noise(1,i)^2+noise(2,i)^2,noise(1,i)^2;noise(1,i)^2,noise(3,i)^2+noise(1,i)^2]);
        q(i+4)=1-mvncdf([0,0],...
        [firing_rate_mean(1,i)-firing_rate_mean(3,i),firing_rate_mean(2,i)-firing_rate_mean(3,i)],...
        [noise(3,i)^2+noise(1,i)^2,noise(3,i)^2;noise(3,i)^2,noise(3,i)^2+noise(2,i)^2]);
    end
    
    answer_relative(V3/20+1)=mean(p./q);
    clear p q;
end
figure
V3temp=0:20:200;
plot(V3temp,answer_relative,'-ro');
hold on;
plot(V3temp,answer_absolute,'-ko');
box off;
legend({'relative model','absolute model'},'Location','NorthWest');
legend boxoff;
xlabel('Distractor value');
ylabel('Efficiency');
% print('Figure2C','-dpng','-r600');
%%
function [firing_rate_mean,noise]=fun(V,String)

if(String=="relative")
    K=100;
    sigmaH=50;
    w=1;
    sigmafixed=8;
    S=0;
    firing_rate_mean=K.*V./(sigmaH+sum(w.*V));
    noise=sqrt(sigmafixed^2+S.*V);
elseif(String=="absolute")
    Ka=0.288;
    sigmafixed=8;
    S=0;
    firing_rate_mean=Ka.*V;
    noise=sqrt(sigmafixed^2+S.*V);
end


end