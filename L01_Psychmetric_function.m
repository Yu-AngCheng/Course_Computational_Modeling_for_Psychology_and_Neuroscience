% Parameter fitting of psychophysical curve
I = [5, 7, 8, 10, 12]'; % Intensity of the light
N = [10, 12, 10, 11, 9]'; % Number of trials
C = [2, 3, 7, 9, 9]'; % Number of positive trials

LB = [0, 0]; % Lower bound
UB = [Inf, Inf]; % Upper bound
x0 = [rand, rand];% Initial seed
[paramsEst, minuslli, exitflag] = ...
fminsearchbnd(@(params)PsychFun(params, I, N, C), x0, LB, UB);
alphaHat = paramsEst(1);
betaHat = paramsEst(2);
lli = - minuslli;


function minuslli = PsychFun(params, I, N, C)
alpha = params(1);
beta = params(2);
% pc is the probability of making positive choices 
pc = 1-exp(-(I/alpha).^beta); % psychometric function;
pc(pc<1e-16) = 1e-16; % to avoid numerical problems near 0 or 1 because of log function
pc(pc>1-1e-16) = 1-1e-16;% to avoid numerical problems near 0 or 1 because of log function
minuslli = -sum(C.*log(pc)+(N-C).*log(1-pc));
end