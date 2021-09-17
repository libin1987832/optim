% The best codes handle N = 20,000 as long as the matrix is very sparse.
% N   = 3000; M = 4000; % Large scale. Things start to get interesting
addpath(genpath(pwd))
N   = 500; M = 1500;     % at this size, some algo take a long time!
% N   = 100; M = 150;     % at this size, all algorithms take < 14 seconds
A   = randn(M,N);
A =[A,-eye(M)];
N = N+M;
b   = randn(M,1);

fcn     = @(x) norm( A*x - b)^2;
% here are two equivalent ways to make the gradient. grad2 is sometimes faster
grad1    = @(x) 2*A'*(A*x-b);
AtA     = A'*A; Ab = A'*b;
grad2    = @(x) 2*( AtA*x - Ab );

grad    = grad2;


x = [];
time = [];
%%
l  = zeros(N,1);    % lower bound
u  = inf(N,1);      % there is no upper bound
tstart=tic;
fun     = @(x)fminunc_wrapper( x, fcn, grad);
% Request very high accuracy for this test:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 5;
if N > 10000
    opts.m  = 50;
end
% Run the algorithm:
[xk, ~, info] = lbfgsb(fun, l, u, opts );
t=toc(tstart)
% Record results
x.lbfgsb    = xk;
time.lbfgsb = t;
% 
% %%
% if exist( 'itron.m', 'file' )
% 
%     x0   = zeros(N,1);
%     xl   = zeros(N,1);
%     xu   = +1e300*ones(N,1);
%     fmin = -1e300;
%     H       = sparse(AtA/2); % will crash if not a sparse matrix
%     tstart=tic;
%     hess    = @(x) H;
%     fun     = @(x)fminunc_wrapper( x, fcn, grad, hess );
%     [xk, fval, exitflag, output] = itron(fun, x0, xl, xu, fmin );
%     t=toc(tstart)
%     x.tron    = xk;
%     time.tron = t;
% end

%%
% if exist( 'activeset.m', 'file' )
%     tstart=tic;
%     [xk,y]  = activeset(A,b);
%     t=toc(tstart)
%     x.activeset    = xk;
%     time.activeset = t;
% end
% %%
% if exist( 'blocknnls.m', 'file' )
%     tstart=tic;
%     [xk]  = blocknnls(A,b, 'fixed');
%     t=toc(tstart)
%     x.blockPivot    = xk;
%     time.blockPivot = t;
% end
%%
if exist( 'newton.m', 'file' ) && N < 500
    tstart=tic;
    [xk,y]  = newton(A,b, ones(N,1), 100); % can't have 0 starting vector
    t=toc(tstart)
    x.newton    = xk;
    time.newton = t;
else
    fprintf('Skipping Newton method because we can''t find it, or it is too slow\n');
end
%%
if exist( 'pcnnls.m', 'file' ) && N < 500
    tstart=tic;
    [xk,y,nits]  = pcnnls(A,b,ones(N,1), 3000);
    t=toc(tstart)
    x.predCorr    = xk;
    time.predCorr = t;
else
    fprintf('Skipping predCorr method because we can''t find it, or it is too slow\n');
end
%%
%
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options = optimoptions('Algorithm','interior-point','TolX',1e-13)
options.Display = 'off';
% options.StepTolerance = 1e-13;
options.OptimalityTolerance = 1e-13;
% options.ConstraintTolerance = 1e-13;
options.MaxIterations = 600;
tstart=tic;
[xk,f1,residual,exitflag,output,ff]=lsqlin(A,b,[],[],[],[],l,u,zeros(N,1),options);
t=toc(tstart)
x.lsqlin   = xk;
time.lsqlin = t; 

%%
% tstart=tic;
% xk = lsqnonneg(A,b);
% t=toc(tstart)
% x.lsqnonneg   = xk;
% time.lsqnonneg = t; 
%%
if exist( 'fnnls.m', 'file' )
    tstart=tic;
    [xk]  = fnnls(A'*A,A'*b);
    t=toc(tstart)
    x.fnnls    = xk;
    time.fnnls = t;
end
%%
if exist( 'solnls.m', 'file' )
 %   opt = solopt;
    opt.maxtime     = 2000;
    opt.verbose     = 0;
    tstart=tic;
    % run their 'BB' variant
    opt.algo = 'BB';
    out     = solnls( A, b, zeros(N,1), opt );
    t=toc(tstart)
    x.PQN_BB   = out.x;
    time.PQN_BB = t;

    % and run their 'PLB' variant (their 'PQN' variant is much slower)
    %   which uses L-BFGS (not to be confused with L-BFGS-B)
    opt.algo = 'PLB';

    tstart=tic;
    out     = solnls( A, b, zeros(N,1), opt );
    t=toc(tstart)

    x.PQN   = out.x;
    time.PQN = t;

end
%%
fMin = Inf;
for f=fieldnames(x)',
    if fcn(x.(f{1})) < fMin,
        fMin = fcn(x.(f{1}));
        best = f{1};
    end
end
xReference = x.(best);
errFcn      = @(x) norm(x-xReference)/norm(xReference);

% Print out info. Verify that the solution is indeed non-negative (hence the
%   min(x) information), and the objective function, and the error
%   against the reference solution. Also display the time.
fprintf('== Size of problem is %d x %d == \n', M, N );
for f=fieldnames(x)',
    fprintf('%10s:  obj is %7.2f, min(x) is %7.1d, err is %.2e, time is %6.3f s\n', ...
        f{1}, fcn(x.(f{1})), min(x.(f{1})), errFcn(x.(f{1})), time.(f{1}) );
end