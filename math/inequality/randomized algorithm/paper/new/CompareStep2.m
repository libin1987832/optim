clear
clc
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

b = 2 * rand(m , 1)-1;
A=real(A);
debug = 1;

tol=1e-5;
maxit_Rand =2000000;
alphaA=[0.5,1,1.5,2];
error=zeros(4,100000);
x_exact=[];
debug=1;
x0=zeros(n,1);
for i = 1:4
    alpha=alphaA(i);
t=clock;
% [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,alpha,maxit_Rand,tol,x_exact,debug);
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = randGuassSeidelNE(A, b, x0,alpha,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b - A * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);
error(i,1:size(error_GS,2))=error_GS;
end
iter=100;
semilogy(error(1,1:iter), 'r') % Matrix A
xlabel('迭代次数') 
ylabel('梯度的误差') 
hold on 
semilogy(error(2,1:iter), 'k') % Matrix M
semilogy(error(3,1:iter), 'g') % Matrix M
semilogy(error(4,1:iter), 'b') % Matrix M
legend('\alpha=0.5', '\alpha=1', '\alpha=1.5', '\alpha=2' )
title('步长对算法的影响(RCD算法求解信号恢复矩阵线性不等式方程组)')

