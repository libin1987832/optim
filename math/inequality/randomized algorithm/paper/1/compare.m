clear
clc
debug = 0;
%% 产生问题矩阵
% 随机矩阵
% m = 500;
% n = 200;

r = 700;
m = 8000;
n = 2*r+1;

%% generate nodes
%tj are drawing randomly from a uniform distribution in [0,1]
t = rand(m,1);

%then ordering by magnitude (since they are all positive, just sort)
t = sort(t);
% just to assign the size
w = zeros(m,1); 
x = zeros(n,1);
%% generate x 
realx = randn(n,1);
imgx = randn(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
%% generate A and b
A = zeros(m,n);
for j = 1:m
    % dealing with special cases when reach the endpoints(nodes)
    if j == 1
        w(j) = (t(2)-t(end)-1)/2;
    elseif j == m
            w(j) = (t(1)+1-t(m-1))/2;
        else
            w(j) = (t(j+1)-t(j-1))/2;
    end              
    for k = -r:r
       A(j,k+r+1) =sqrt(w(j))*exp(2*pi*1i*k*t(j));
    end
end  

%  A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% A=[A;-A];
% b=[b;-b];
% m=2*m;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];
tol=1e-1;

fprintf('%s & %d & %d \n','矩阵维数', m, n);
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



%% GuassSeidel
maxit_Rand =2000000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,1.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'simpleGuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'randGuassSeidel', r_GS, g_GS, iter_GS, tf_GS);



%% 参数的设定
 maxit_IFM = 200;


t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A, b, x0, maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b - A * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

%% FM
% maxit =100;

alpha=1;
maxit_gs=n;
t=clock;
[x_FM,iter_FM,error_k,iter_dFM,index_k] = DFM(A, b, x0, maxit_IFM,alpha,maxit_gs,tol, x_exact,debug);
tf_FM=etime(clock,t);
r = b - A * x_FM;
r(r<0) = 0;
r_FM = norm(r);
g_FM = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM,iter_FM,tf_FM);

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
% maxIter = 1;
% t=clock;
% [x_GS,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
% tf_GS=etime(clock,t);
% r = b - A * x_GS;
% r(r<0) = 0;
% r_GS = norm(r);
% g_GS = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n', 'han', r_GS, g_GS, countFMh, tf_GS);
