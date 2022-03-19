clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
% m = 500;
% n = 200;

r = 50;
m = 700;
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
%clear
% load('dd')
%A=real(A);
debug = 1;
 A=[A;-A];
 realx = randn(n,1);
imgx = randn(n,1);
x=zeros(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
b=A*x;
% b=[b;-b];
 m=2*m;
 A=2*rand(m,n)-1;
 b=2*rand(m,1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];
tol=1e-5;

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
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_SGS,iter_GS,error_SGS,xA_GS,index_GS] = simpleGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_SGS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'simpleGuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_GS,error_RGS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'randGuassSeidel', r_GS, g_GS, iter_GS, tf_GS);


semilogy(error_GS, 'r') % Matrix A
xlabel('迭代次数') 
ylabel('梯度的误差') 
hold on 
semilogy(error_SGS, 'k') % Matrix M
 semilogy(error_RGS, 'g') % Matrix M
legend('CCD', 'UCD', 'RCD' )
title('三种指标选择模式对算法的影响')



%%%%%%%%%%%%
% 矩阵维数 & 8000 & 1401 
% 最开始的目标函数和梯度 & 36.9619 & 28.2101 
% IFM解的目标函数值和梯度   & 7.83054 & 1.8207 
% & GuassSeidel & 6.85029 & 2.84068 & 43431 & 4.646 \\
% & simpleGuassSeidel & 4.65163 & 0.883448 & 29421 & 3.216 \\
% & randGuassSeidel & 4.64761 & 0.883935 & 25218 & 2.868 \\
% & IFM & 4.02986 & 0.623185 & 29 & 5.687 \\
% & FM & 9.11728 & 4.12671 & 33 & 4.011 \\

% 矩阵维数 & 700 & 101 
% 最开始的目标函数和梯度 & 10.7863 & 6.7381 
% IFM解的目标函数值和梯度   & 2.99591 & 0.633142 
% & GuassSeidel & 0.209602 & 0.00789479 & 28886 & 0.131 \\
% & simpleGuassSeidel & 0.368795 & 0.0146253 & 16766 & 0.156 \\
% & randGuassSeidel & 0.258481 & 0.00952142 & 19897 & 0.164 \\
