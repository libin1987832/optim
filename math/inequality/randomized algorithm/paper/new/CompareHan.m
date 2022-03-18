clear
clc
debug = 1;
%% 产生问题矩阵
m = 100;
n = 10;
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
[x_Han,r_Han,countFM_Han,error_Han,beginNW_Han,tf_Han,vk_Han,rkArr_Han]=han(x0,A,b,maxit_Han,x_exact,debug);
tf_PC=etime(clock,t);
r = b - A * x_Han;
r(r<0) = 0;
r_PC = norm(r);
g_PC = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','Han', r_PC, g_PC,countFM_Han,tf_PC);

%% GuassSeidel
maxit_Rand =500;
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
[x_GS,iter_SGS,error_SGS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'UCD', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_RGS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RCD', r_GS, g_GS, iter_GS, tf_GS);
% clear
% clc
% load("data_han")
% figure (1)
% semilogy(1:countFM_Han,error_Han(1:countFM_Han),'k--');
% title('比较坐标下降算法和Han算法的性能')
% ylabel('梯度的范数') 
% xlabel('迭代次数（对于坐标下降法，每n次迭代算一次迭代）') 
% hold on
% iters=250;
% semilogy(1:iters,error_GS(1:iters),'b--');
% semilogy(1:iters,error_SGS(1:iters),'r--');
% semilogy(1:iters,error_RGS(1:iters),'g--');
% 
% % semilogy(1:size(error_GS,2),error_GS,'b--');
% % semilogy(1:size(error_SGS,2),error_SGS,'r--');
% % semilogy(1:size(error_RGS,2),error_RGS,'g--');
% hold off
% legend('Han','CCD','UCD','RCD')

% 矩阵维数 & 5000 & 1000 tol=-1
% 最开始的目标函数和梯度 & 28.9627 & 523.841 
% & Han & 22.3306 & 7.61694e-13 & 7 & 52.975 \\
% & CCD & 22.3306 & 1.35664e-10 & 1000000 & 30.099 \\
% & UCD & 22.3306 & 9.47045e-11 & 1000000 & 31.492 \\
% & RCD & 22.3306 & 9.28481e-11 & 1000000 & 35.044 \\