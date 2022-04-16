clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
% m = 500;
% n = 200;


m = 100;
n = 20;
  A = 2 * rand(m , n)-1;
   b = 2 * rand(m , 1)-1;

x0 = zeros(n , 1);

tol=1e-5;

fprintf('%s & %d & %d \n','矩阵维数', m, n);
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
%% GuassSeidel
maxit_Rand =8000000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CCD', r_GS, g_GS, iter_GS, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_GS,iter_GS,error_SGS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'UCD', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_GS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RCD', r_GS, g_GS, iter_GS, tf_GS);


% semilogy(error_GS, 'r') % Matrix A
xlabel('Iterations') 
ylabel('Error') 
hold on 
semilogy(error_GS, 'k') % Matrix M
semilogy(error_SGS, 'g') % Matrix M
semilogy(error_RGS, 'b') % Matrix M
legend('CCD', 'UCD', 'RCD' )

Pi=eye(m)-