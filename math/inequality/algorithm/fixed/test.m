addpath('../../dataInequality')
addpath('../../algorithmInequality')
clear
clc
m = 1000;
n = 500;
[A,b,x00] = randInequality(m,n,2,-1);
%tic
[Q,R] = qr(A);
x0=x00;
rk=b-A*x0; 
[xk,rk] = FMQR(x0,Q,R,A,b,rk);
iter = 120;
for i = 1:iter
    [xk,rk] = FMQR(x0,Q,R,A,b,rk);
    x0=xk;
end
%toc
[rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
nRkN
nGARk
%tic
x0=x00;
iter = 150;
for i = 1:iter
    [xk,rk]=FMLSQR(x0,A,b,rk);
    x0=xk;
end
%toc
[rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
nRkN
nGARk