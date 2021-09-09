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