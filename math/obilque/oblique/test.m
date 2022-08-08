clc;clear;close all;
%³õÊ¼Öµ
m=300;n=50;
A=rand(m,n);
b=rand(m,1);RRE=0.5e-6;x_0=zeros(n,1);

if m>=n
    x_exact=(A'*A)\A'*b;
    b_A=b-A*x_exact;
else
    b_A=b-A'\(A'*b);
end
%%
debug=1;
tol=1e-10;
iter=1e5;
[xC,errorC]=GSO(A,b,RRE,x_0,tol,iter,b_A,debug);
[xR,errorR]=RGSO2(A,b,RRE,x_0,tol,iter,b_A,debug);
[xGS,errorGS] = GS(A, b, RRE,x_0,1,iter,b_A,debug);
iters=1000;
semilogy(1:iters,errorC(1:iters),'b.');
hold on
semilogy(1:iters,errorR(1:iters),'r.');
semilogy(1:iters,errorGS(1:iters),'g.');



