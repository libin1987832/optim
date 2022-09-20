clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵
m = 10000;
n = 3000;
A = 2*rand(m,n)-1;
b = rand(m,1);

x0 = zeros(n , 1);



x_exact=[];
%% 参数的设定
 maxit_IFM = 150;
maxit_LSQR = 3;
tol=1e-1;
% tol=[];

t=clock;
[x_IFMs,iter_IFMs,error_IFMs,xA_IFMs,index_IFMs] = IFMs(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFMs=etime(clock,t);
r = b - A * x_IFMs;
r(r<0) = 0;
r_IFMs = norm(r);
g_IFMs = norm(A'*r);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFMs, g_IFMs, iter_IFMs*(m+n)*maxit_LSQR,tf_IFMs);
t=clock;
[x_IFMp,iter_IFMp,error_IFMp,xA_IFMp,index_IFMp] = IFMp(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFMp=etime(clock,t);
r = b - A * x_IFMp;
r(r<0) = 0;
r_IFMp = norm(r);
g_IFMp = norm(A'*r);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFMp, g_IFMp,iter_IFMp*(m+n)*maxit_LSQR,tf_IFMp);

