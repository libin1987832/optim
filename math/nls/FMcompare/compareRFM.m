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
 maxit_FM = 1000;
maxit_gs = 1;
tol=1e-1;
% tol=[];
alpha = 1;
t=clock;
[x_FM,iter_FM,error_k,iter_k,index_k] = DFM(A, b, x0, maxit_FM,alpha,maxit_gs,tol, x_exact,debug);
tf_FM = etime(clock,t);
r = b - A * x_FM;
r(r<0) = 0;
r_FM = norm(r);
g_FM = norm(A'*r);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM, iter_FM*(m+n)*maxit_gs, tf_FM);


