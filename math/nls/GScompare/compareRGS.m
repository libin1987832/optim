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
maxit_GS =70000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,1.0,maxit_GS,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

