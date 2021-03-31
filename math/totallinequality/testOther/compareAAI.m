addpath('../debug')

[A,b,x0] = readData(1,400,400,200);

[m,n] =size(A);
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options = optimoptions('Algorithm','interior-point','TolX',1e-13)
options.Display = 'off';
% options.StepTolerance = 1e-13;
options.OptimalityTolerance = 1e-15;
options.ConstraintTolerance = 1e-20;
options.MaxIterations = 200;


% load('testMIL')
% options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
[m,n] =size(A);

tic;
[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
    [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);
toc;
tic;
[x1,f1,residual,exitflag,output,ff] = lsqlin(A,b,...
    [],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
toc;

