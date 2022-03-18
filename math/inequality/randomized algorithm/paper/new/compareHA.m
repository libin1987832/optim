clear
clc
debug = 0;
%% 产生问题矩阵
m = 2000;
n = 200;
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;

x0 = zeros(n , 1);
tol=1e-50-1;
fprintf('%s & %d & %d tol=%g\n','矩阵维数', m, n, tol);
%% 基于IFM的算法找到一个解
maxit_LSQR = 3;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);
% [x_exact, ~, ~, ~, ~] = IFM(A, b, x0,10, maxit_LSQR , 1e-10,[],debug);
% r = b - A * x_exact;
% r(r<0) = 0;
% norm_rexact = norm(r);
% norm_gexact = norm(A'*r);
% fprintf('%s & %g & %g \n','IFM解的目标函数值和梯度  ', norm_rexact, norm_gexact);
 x_exact=[];
%% 参数的设定
maxit_Han = 50;
z0=zeros(m,1);
t=clock;
[x_Han,r_Han,countFM_Han,error_Han,beginNW_Han,tf_Han,vk_Han,rkArr_Han]=Hybrid(x0,A,b,maxit_Han,1e-11,x_exact,debug);
tf_PC=etime(clock,t);
r = b - A * x_Han;
r(r<0) = 0;
r_PC = norm(r);
g_PC = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','Han', r_PC, g_PC,countFM_Han,tf_PC);

%% GuassSeidel
maxit_Rand =5000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand*2,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CCD', r_GS, g_GS, iter_GS, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_GS,iter_SGS,error_SGS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand*5,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'UCD', r_GS, g_GS, iter_SGS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_RGS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand*5,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RCD', r_GS, g_GS, iter_RGS, tf_GS);
