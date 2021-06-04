addpath('../dataInequality/');
clc
clear
m = 1000;
n = 500;
rangeMax = 2;
rangeMin = -2;
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x00 = zeros(n , 1);
x0=x00;
rk0 = b-A*x0;
iter = 3000;
[rk, rkN, nRkN, nGARk] = residual(A,b,x0);
fprintf('residual %g gradient %g \n',nRkN,nGARk);
[nf,iffind,xk,sumrkA]=generateNf(A,b,x0,iter);
%[nf,iffind]=generateNf(A,b,x00,100);
[rk, rkN, nRkN, nGARk] = residual(A,b,xk);
fprintf('residual %g gradient %g ',nRkN,nGARk);
sumrkA(1,iter-100:iter)


% rkA=zeros(m,iter);
% sumrkA=zeros(1,iter);
% rkA(:,1)=rk0;
% sumrkA=sum(rk0>0);
% for i = 2:iter
% [xk,rk]=Lei(x0,A,b,3,rk0);
% x0=xk;
% rk0=rk;
% rkA(:,i)=rk;
% sumrkA(:,i)=sum(rk>0);
% end
% sumrkA