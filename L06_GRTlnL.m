function [lnL,predP]=GRTlnL(theta,x,data,N)
bound1=theta(1);
bound2=theta(2);
sd=theta(3);
predP=normcdf((bound1-x)./sd)+1-normcdf((bound2-x)./sd);
lnL=-sum(data.*log(predP)+(N-data).*log(1-predP));