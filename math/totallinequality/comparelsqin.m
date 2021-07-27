addpath('./IPG')
addpath('./exact')
addpath('./hybrid')
addpath('./subproblem')
addpath('./GNP')
addpath('./exponent')
addpath('./hybridfast')
addpath('../inequality')
addpath('../inequality/algorithmInequality')
[A,b,x0] = readData(2,1000,1000,100);
[m,n] =size(A);
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options = optimoptions('Algorithm','interior-point','TolX',1e-13)
options.Display = 'off';
% options.StepTolerance = 1e-13;
options.OptimalityTolerance = 1e-13;
% options.ConstraintTolerance = 1e-13;
options.MaxIterations = 600;
tic;
[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
    [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);
tf0=toc;
xk0 = x1(1:n);
[rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, xk0 , [], 1); 
fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ls',normr0,xmin0,normKKT0,min(Ar0),tf0);