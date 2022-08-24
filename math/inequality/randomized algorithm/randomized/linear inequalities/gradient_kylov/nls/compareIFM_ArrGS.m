clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
m = 3000;
n = 300;

A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];


%% 基于IFM的算法找到一个解
maxit_LSQR = 3;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);
[x_exact, ~, ~, ~, ~] = IFM(A, b, x0,10, maxit_LSQR , 1e-10,[],debug);
r = b - A * x_exact;
r(r<0) = 0;
norm_rexact = norm(r);
norm_gexact = norm(A'*r);
fprintf('%s & %g & %g \n','IFM解的目标函数值和梯度  ', norm_rexact, norm_gexact);
x_exact=[];
%% 参数的设定
 maxit_IFM = 150;
% 
tol=1e-1;
% tol=[];

t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM*(m+n)*maxit_LSQR,tf_IFM);


%% GuassSeidel
maxit_GS =35000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,1.0,maxit_GS,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

maxit_FM = 1000;
alpha = 1;
maxit_gs = 1;
t=clock;
[x_FM,iter_FM,error_k,iter_k,index_k] = DFM(A, b, x0, maxit_FM,alpha,maxit_gs,tol, x_exact,debug);
tf_FM=etime(clock,t);
r = b - A * x_FM;
r(r<0) = 0;
r_FM = norm(r);
g_FM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'FM', r_FM , g_FM, iter_FM*n, tf_FM);

%% Rand
maxit_Rand =35000;
t=clock;
[x_RGS,iter_RGS,error_RGS,xA_RGS,index_RGS] = randGSNE(A, b, x0,1.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_RGS;
r(r<0) = 0;
r_RGS = norm(r);
g_RGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RGS', r_GS, g_GS, iter_GS, tf_GS);
