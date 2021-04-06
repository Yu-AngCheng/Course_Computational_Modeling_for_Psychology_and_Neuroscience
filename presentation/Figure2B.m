clear
%% absolute
set(0,'defaultfigurecolor','w') 
h=figure;
v=VideoWriter('Figure2B_absolute','MPEG-4');
v.FrameRate=1;
v.Quality=100;
open(v)
for V3temp=0:20:200
    V1=[110,120,130,140,150,160,170,180,190];
    V2=[150,150,150,150,150,150,150,150,150];
    V3=repmat(V3temp,size(V1));
    V=[V1;V2;V3];
    [firing_rate_mean,noise]=fun(V,"absolute");
    for i=1:length(V1)
        p(i)=mvncdf([0,0],...
            [firing_rate_mean(2,i)-firing_rate_mean(1,i),firing_rate_mean(3,i)-firing_rate_mean(1,i)],...
            [noise(1,i)^2+noise(2,i)^2,noise(1,i)^2;noise(1,i)^2,noise(3,i)^2+noise(1,i)^2]);
        q(i)=1-mvncdf([0,0],...
            [firing_rate_mean(1,i)-firing_rate_mean(3,i),firing_rate_mean(2,i)-firing_rate_mean(3,i)],...
            [noise(3,i)^2+noise(1,i)^2,noise(3,i)^2;noise(3,i)^2,noise(3,i)^2+noise(2,i)^2]);
    end
    plot(V1-V2,p./q);
    legend({'V3=0','V3=20','V3=40','V3=60','V3=80','V3=100',...
        'V3=120','V3=140','V3=160','V3=180','V3=200'});
    legend boxoff;
    legend('Location','northwest');
    ylim([0,1]);
    box off
    xlabel('Target value difference');
    ylabel('Relative target choice');
    title('Absolute model');
    writeVideo(v,getframe(gcf))
    %     Fig1(V3temp/20+1)=getframe(gca);
    hold on;
    clear p q;
end
close(v);
% print('Figure2B_absolute','-dpng','-r600');
%% normalization
h=figure;
v=VideoWriter('Figure2B_relative','MPEG-4');
v.FrameRate=1;
v.Quality=100;
open(v)
for V3temp=0:20:200
    V1=[110,120,130,140,150,160,170,180,190];
    V2=[150,150,150,150,150,150,150,150,150];
    V3=repmat(V3temp,size(V1));
    V=[V1;V2;V3];
    [firing_rate_mean,noise]=fun(V,"relative");
    for i=1:length(V1)
        p(i)=mvncdf([0,0],...
            [firing_rate_mean(2,i)-firing_rate_mean(1,i),firing_rate_mean(3,i)-firing_rate_mean(1,i)],...
            [noise(1,i)^2+noise(2,i)^2,noise(1,i)^2;noise(1,i)^2,noise(3,i)^2+noise(1,i)^2]);
        q(i)=1-mvncdf([0,0],...
            [firing_rate_mean(1,i)-firing_rate_mean(3,i),firing_rate_mean(2,i)-firing_rate_mean(3,i)],...
            [noise(3,i)^2+noise(1,i)^2,noise(3,i)^2;noise(3,i)^2,noise(3,i)^2+noise(2,i)^2]);
    end
    plot(V1-V2,p./q);
    legend({'V3=0','V3=20','V3=40','V3=60','V3=80','V3=100',...
        'V3=120','V3=140','V3=160','V3=180','V3=200'});
    legend boxoff;
    legend('Location','northwest');
    ylim([0,1]);
    box off
    xlabel('Target value difference');
    ylabel('Relative target choice');
    title('Relative model');
    writeVideo(v,getframe(gcf))
    hold on;
    clear p q;
end
close(v)
% print('Figure2B_relative','-dpng','-r600');
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