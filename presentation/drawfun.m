function drawfun(handle,miu,sigmaFix,S,color)
    sigma=sqrt(sigmaFix^2+S*miu);
    x=linspace(miu-4*sigma,miu+4*sigma,5000);
    y=normpdf(x,miu,sigma);
    plot(handle,x,y,color);
end