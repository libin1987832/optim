clear
clc
m=5000;
n=500;
A=2*rand(m,n)-1;
b=2*rand(m,1)-1;
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
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A, b, x0,alpha,maxit_Rand,tol,x_exact,debug);
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
title('步长对算法的影响(CCD算法求解均匀矩阵构造线性不等式方程组)')

