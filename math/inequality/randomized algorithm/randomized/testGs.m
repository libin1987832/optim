clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
m = 1000;
n = 200;

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
% maxit_IFM =210;
% 
tol=1e-10;
tol=[];
% 
% t=clock;
% [x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
% tf_IFM=etime(clock,t);
% r = b - A * x_IFM;
% r(r<0) = 0;
% r_IFM = norm(r);
% g_IFM = norm(A'*r);
% fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
% fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);


%% GuassSeidel
maxit_Rand =1000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,1,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

%% GaussSeidel
% [U,S,V]=eig(A);
%%%
maxit_Rand =1000;
t=clock;
% [x_WGS,iter_WGS,error_WGS,xA_WGS,index_WGS] = wrandomizedGaussSeidelNE(A, b, x0,20,10, maxit_Rand, tol,x_exact,debug);
 [x_WGS,iter_WGS,error_WGS,xA_WGS,index_WGS] = GuassSeidelNE(A, b, x0,0.8,maxit_Rand,tol,x_exact,debug);
% [x_WGS,iter_WGS,error_WGS,xA_WGS,index_WGS] = randomizedGaussSeidelNE(A, b, x0,2,maxit_Rand, tol,x_exact,debug);

tf_WGS=etime(clock,t);
r = b - A * x_WGS;
r(r<0) = 0;
r_WGS = norm(r);
g_WGS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'weight Gauss', r_WGS, g_WGS, iter_WGS, tf_WGS);
%% 画图
if debug
figure
% h=semilogy(xA_IFM, error_IFM, 'k.');
% h.LineStyle = '--';
display=1:1:maxit_Rand;
h=semilogy(xA_GS(display), error_GS(display), 'r+');
h.LineStyle = '--';
hold on
h=semilogy(xA_WGS(display), error_WGS(display), 'b*');
h.LineStyle = '--';
legend('Gauss Seidel','Rand Guass');
xlabel('the iterative numbers');
ylabel('the norm of the gradient');
end

