V3=0:1:160;
V1=repmat(150,1,length(V3));
V2=repmat(140,1,length(V3));
V=[V1;V2;V3];
iteration=10000;
%% absolute
firing_rate_with_error=fun(V,"absolute",iteration);
[~,Index]=max(firing_rate_with_error);
Index=squeeze(Index)';
p1_over_p2_arbitary=sum(Index==1)./sum(Index==2);
semilogy(V3,p1_over_p2_arbitary,'ko')
hold on;
%% normalization
firing_rate_with_error=fun(V,"relative",iteration);
[~,Index]=max(firing_rate_with_error);
Index=squeeze(Index)';
p1_over_p2_normalization=sum(Index==1)./sum(Index==2);
semilogy(V3,p1_over_p2_normalization,'ko','MarkerFaceColor','k');
ylim([10 1000])
%%
legend({'Absolute','Relative'},'Location','Northwest');
xlabel('Distractor value');
ylabel('Target choice ratio');
print('Figure1C','-dpng','-r600');
%%
function firing_rate_with_error=fun(V,String,iteration)

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
firint_rate_mean=repmat(firing_rate_mean,1,1,iteration);
noise=normrnd(0,sigmafixed,size(firint_rate_mean))+normrnd(0,repmat(sqrt(S.*V),1,1,iteration),size(firint_rate_mean));
firing_rate_with_error=firing_rate_mean+noise;

end