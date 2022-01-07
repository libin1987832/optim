clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵
 m = 5000;
 n = 2000;
randone= rand(1,n)<0.5;

 A = 2 * rand(m , n)-1;
A(:,randone)=A(:,randone)/1000;
b = 2 * rand(m , 1)-1;
A=[A;-A];
b=[b;-b];
m=2*m;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];
tol=1e-1;


%% 基于IFM的算法找到一个解
maxit_LSQR = 3;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);
x_exact = [];


%% GuassSeidel
maxit_Rand =2000000;
t=clock;
 [x_CGS,iter_CGS,~,~,~] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_CGS=etime(clock,t);
r = b - A * x_CGS;
r(r<0) = 0;
r_CGS = norm(r);
g_CGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CGS($\lambda=2$)', r_CGS, g_CGS, iter_CGS, tf_CGS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_SGS,iter_SGS,~,~,~] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_SGS=etime(clock,t);
r = b - A * x_SGS;
r(r<0) = 0;
r_SGS = norm(r);
g_SGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'SGS($\lambda=2$)', r_SGS, g_SGS, iter_SGS, tf_SGS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_RGS,iter_RGS,~,~,~] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_RGS=etime(clock,t);
r = b - A * x_RGS;
r(r<0) = 0;
r_RGS = norm(r);
g_RGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'RGS($\lambda=2$)', r_RGS, g_RGS, iter_RGS, tf_RGS);



%% 参数的设定
 maxit_IFM = 200;


t=clock;
[x_IFM,iter_IFM,~,~,~] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% FM

alpha=1;
maxit_gs=n;
t=clock;
[x_FM,iter_FM,~,~,~] = DFM(A, b, x0, maxit_IFM,alpha,maxit_gs,tol, x_exact,debug);
tf_FM=etime(clock,t);
r = b - A * x_FM;
r(r<0) = 0;
r_FM = norm(r);
g_FM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM,iter_FM,tf_FM);

