clear
clc
debug = 1;
%% 产生问题矩阵
% 随机矩阵
% m = 1000;
% n = 200;

r = 500;
m = 10000;
n = 2*r+1;

% generate nodes
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

% A = 2 * rand(m , n)-1;
A = [A;-A];
b = 2 * rand(m , 1)-1;
b = [b;-b];
% b=A*ones(n,1);
x0 = zeros(n , 1);
% save('test2.mat','A','b','x0')
%load('test.mat');
% 二维矩阵
% A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
% 不一致情况下的正解
% x_exact=[1/2;1/3];
tol=[];


%% 基于IFM的算法找到一个解
maxit_LSQR = 3;
r = b - A * x0;
r(r<0) = 0;
norm_r0 = norm(r);
norm_g0 = norm(A'*r);
fprintf('%s & %g & %g \n','最开始的目标函数和梯度', norm_r0, norm_g0);

x_exact=[];

tol = 1e-1;

%% GuassSeidel
maxit_Rand =50000;
t=clock;
 [x_GS,iter_GS1,error_GS1,xA_GS1,index_GS] = GuassSeidelNE(A, b, x0,1.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CGS(\lambda=1)', r_GS, g_GS, iter_GS1, tf_GS);


%% simpleGuassSeidel
maxit_Rand =50000;
t=clock;
 [x_GS,iter_GS2,error_GS2,xA_GS2,index_GS] = GuassSeidelNE(A, b, x0,1.0+min(m/n,n/m),maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CGS(\lambda=1+min(m/n,n/m))', r_GS, g_GS, iter_GS2, tf_GS);

%% randGuassSeidel
maxit_Rand =50000;
t=clock;
 [x_GS,iter_GS3,error_GS3,xA_GS3,index_GS] = GuassSeidelNE(A, b, x0,1.5,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CGS(\lambda=1.5)', r_GS, g_GS, iter_GS3, tf_GS);

maxit_Rand =50000;
t=clock;
 [x_GS,iter_GS4,error_GS4,xA_GS4,index_GS] = GuassSeidelNE(A, b, x0,2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'CGS(\lambda=2)', r_GS, g_GS, iter_GS4, tf_GS);

%% 画图
if debug
figure
% h=semilogy(xA_IFM, error_IFM, 'k.');
% h.LineStyle = '--';
 sum = min([iter_GS1,iter_GS2,iter_GS3,iter_GS4]) ;
iter_GS1 = sum;
iter_GS2 = sum;
iter_GS3 = sum;
iter_GS4 = sum;
display=1:1:iter_GS1;
h=semilogy(xA_GS1(display), error_GS1(display), 'r+');
h.LineStyle = '--';
hold on
display=1:1:iter_GS2;
h=semilogy(xA_GS2(display), error_GS2(display), 'b+');
h.LineStyle = '--';
display=1:1:iter_GS3;
h=semilogy(xA_GS3(display), error_GS3(display), 'g+');
h.LineStyle = '--';
display=1:1:iter_GS4;
h=semilogy(xA_GS4(display), error_GS4(display), 'k+');
h.LineStyle = '--';
legend('\lambda = 1','\lambda = 1+min(m/n,m/m)','\lambda = 1.5','\lambda = 2');
xlabel('the iterative numbers');
ylabel('the norm of the gradient');
end



