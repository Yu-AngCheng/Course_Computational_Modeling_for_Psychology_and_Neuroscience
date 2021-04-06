V3=100;
V1=150;
V2=140;
V=[V1;V2;V3];
%% absolute
[firing_rate_mean,noise]=fun(V,"absolute");
X=0:0.1:50;
P=normpdf(X,firing_rate_mean,noise);
figure;
plot(X,P(1,:),'b');
hold on;
plot(X,P(2,:),'k');
plot(X,P(3,:),'r');
%% normalization
[firing_rate_mean,noise]=fun(V,"relative");
X=0:0.1:50;
P=normpdf(X,firing_rate_mean,noise);
figure;
plot(X,P(1,:),'b');
hold on;
plot(X,P(2,:),'k');
plot(X,P(3,:),'r');
%%
function [firing_rate_mean,noise]=fun(V,String)

if(String=="relative")
    K=100;
    sigmaH=50;
    w=1;
    sigmafixed=1;
    S=0;
elseif(String=="absolute")
    K=100/(50+140+150);
    sigmaH=1;
    w=0;
    sigmafixed=1;
    S=0;
end
firing_rate_mean=K.*V./(sigmaH+sum(w.*V));
noise=sqrt(sigmafixed^2+S.*V);

end