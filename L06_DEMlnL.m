function [lnL,predP]=DEMlnL(theta,x,feedback,data,N)
c=theta(1);
gamma=theta(2);
for i=1:length(x)
    s=exp(-c.*abs(x(i)-x));
    sumA=sum(s.*feedback(1,:));
    sumB=sum(s.*feedback(2,:));
    predP(i)=(sumA.^gamma)/(sumA.^gamma+sumB.^gamma);
end
lnL=-sum(data.*log(predP)+(N-data).*log(1-predP));