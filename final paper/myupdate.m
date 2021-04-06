function answer=myupdate(IL,RL,AL,CL,ML,string)
if string=="RL"
    lambda=0.95;
    A=...
        [1.5,-1,-1,0,0;
        -1,1.5,-1,0,0;
        -0.25,-0.25,1.25,-0.25,-0.25;
        0,0,-1,1.5,-1;
        0,0,-1,-1,1.5];
    B=...
        [2.0,0.33,0.33,0,0;
        0.33,2.0,0.33,0,0;
        0,0,1.5,0,0;
        0,0,0.33,2.0,0.33;
        0,0,0.33,0.33,2.0];
    ksai=normrnd(0,[0.2,0.2,1.25,0.2,0.2]);
    thetax=0.5;
    X=RL*A+IL*B-thetax+ksai;
    answer=lambda*RL+(1-lambda)*logistic(X,AL);
elseif string=="CL"
    lambda=0.75;
    A=...
        [1,-1;
        -1,1];
    B=...
        [1,0;
        0,1;
        0,0;
        1,0;
        0,1];
    ksai=normrnd(0,[0.025,0.025]);
    thetax=0.85;
    X=CL*A+RL*B-thetax+ksai;
    answer=lambda*CL+(1-lambda)*logistic(X,AL);
elseif string=="ML"
    lambda=0.925;
    A=...
        [0.9,-1;
        -1,0.9];
    B=...
        [1.5,0;
        0,1.5;
        0,0;
        1.5,0;
        0,1.5];
    ksai=normrnd(0,[0.25,0.25]);
    thetax=2;
    X=ML*A+RL*B-thetax+ksai;
    answer=lambda*ML+(1-lambda)*logistic(X,AL);
elseif string=="AL"
    lambdax=0.92;
    lambday=0.996;
    lambdag=0.98;
    E=CL(1)*CL(2);
    input=sum(RL.*[1,1,0,1,1]);% or maybe mean
%     if (E>=0.5)&&((RL(4)>0.62)||(RL(5)>0.62))
    if E>=0.5
        c=1;
    else
        c=3;
    end
    X=lambdax*AL(1)+(1-lambdax)*logistic(c*(2*AL(1)-4*AL(2)+input-1.25),AL);
    Y=lambday*AL(2)+(1-lambday)*logistic(c*(3*AL(1)-1.5),AL);
    G=lambdag*AL(3)+(1-lambdag)*AL(1);
    answer=[X,Y,G];
end
end
