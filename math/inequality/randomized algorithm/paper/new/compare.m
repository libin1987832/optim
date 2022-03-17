clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵
% m = 500;
% n = 200;


m = 8000;
n = 5000;
  A = 2 * rand(m , n)-1;
   b = 2 * rand(m , 1)-1;

x0 = zeros(n , 1);

tol=1e-10;

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
%% 参数的设定
 maxit_IFM = 500;
t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% GuassSeidel
maxit_Rand =3000000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,1.0,maxit_Rand,tol,x_exact,debug);
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








%% FM
% maxit =100;

% alpha=1;
% maxit_gs=n;
% t=clock;
% [x_FM,iter_FM,error_k,iter_dFM,index_k] = DFM(A, b, x0, maxit_IFM,alpha,maxit_gs,tol/50, x_exact,debug);
% tf_FM=etime(clock,t);
% r = b - A * x_FM;
% r(r<0) = 0;
% r_FM = norm(r);
% g_FM = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM,iter_FM,tf_FM);

%% RFM
% maxit =18;
% tol=[];
% alpha=1;
% maxit_R=n;
% t=clock;
% [x_FM,iter_FM,error_k,iter_dFM,index_k] = RFM(A, b, x0, maxit,alpha,maxit_R, tol, x_exact,debug);
% tf_FM=etime(clock,t);
% r = b - A * x_FM;
% r(r<0) = 0;
% r_FM = norm(r);
% g_FM = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n','RFM', r_FM, g_FM,iter_FM,tf_FM);


% %% han
% maxIter = 5;
% t=clock;
% [x_GS,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
% tf_GS=etime(clock,t);
% r = b - A * x_GS;
% r(r<0) = 0;
% r_GS = norm(r);
% g_GS = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n', 'han', r_GS, g_GS, countFMh, tf_GS);


% semilogy(error_GS, 'r') % Matrix A
% xlabel('Iterations') 
% ylabel('Error') 
% hold on 
% semilogy(error_SGS, 'k') % Matrix M
% semilogy(error_RGS, 'g') % Matrix M
% semilogy(error_IFM, 'b') % Matrix M
% legend('CCD', 'UCD', 'RCD','IFM' )

%%%%%%%%%%%%
% 矩阵维数 & 8000 & 1401 
% 最开始的目标函数和梯度 & 36.9619 & 28.2101 
% IFM解的目标函数值和梯度   & 7.83054 & 1.8207 
% & GuassSeidel & 6.85029 & 2.84068 & 43431 & 4.646 \\
% & simpleGuassSeidel & 4.65163 & 0.883448 & 29421 & 3.216 \\
% & randGuassSeidel & 4.64761 & 0.883935 & 25218 & 2.868 \\
% & IFM & 4.02986 & 0.623185 & 29 & 5.687 \\
% & FM & 9.11728 & 4.12671 & 33 & 4.011 \\
