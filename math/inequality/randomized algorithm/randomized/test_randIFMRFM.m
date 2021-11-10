clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
m = 500000;
n = 1000;

A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test.mat','A','b','x0')
load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];


%% 基于IFM的算法找到一个解
maxit_LSQR = 10;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);
[x_exact, ~, ~, ~, ~] = IFM(A, b, x0,1000, maxit_LSQR , 1e-15,[],debug);
r = b - A * x_exact;
r(r<0) = 0;
norm_rexact = norm(r);
norm_gexact = norm(A'*r);
fprintf('%s & %g & %g \n','IFM解的目标函数值和梯度  ', norm_rexact, norm_gexact);
if debug == 0
x_exact=[];
end
%% 参数的设定
maxit_IFM = 100;
maxit_Rand = 100;
tol=1e-5;
tol=[];
%% IFM算法求解问题
t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %s & %s & %s & %s \\\\\n', 'alg', 'norm(r_+)', 'norm(Ar_+)', 'iteration', 'time');
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% RFM算法求解问题
t=clock;
maxit_R =300;
[x_RFM,iter_RFM,error_RFM,xA_RFM,index_RFM] = RFM(A, b, x0, maxit_IFM, maxit_R ,tol, x_exact,debug);
tf_RFM=etime(clock,t);
r = b - A * x_RFM;
r(r<0) = 0;
r_RFM = norm(r);
g_RFM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','RFM', r_RFM, g_RFM,iter_RFM,tf_RFM);


%% 画图
if debug
figure
h=semilogy(1:(iter_IFM+1), error_IFM, 'k.');
h.LineStyle = '--';
hold on
h=semilogy(1:(iter_RFM+1), error_RFM, 'r+');
h.LineStyle = '--';
legend('IFM','RFM');
xlabel('the iterative numbers');
ylabel('the norm of the gradient');
end

