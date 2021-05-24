addpath('../dataInequality/');
addpath('../FM/')
clc
clear
m = 1000;
n = 100;
rangeMax = 2;
rangeMin = -2;
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
rk0 = b-A*x0;
iter = 50;
rkA=zeros(m,iter);
sumrkA=zeros(1,iter);
rkA(:,1)=rk0;
sumrkA=sum(rk0>0);
for i = 2:iter
[xk,rk]=Lei(x0,A,b,3,rk0);
x0=xk;
rk0=rk;
rkA(:,i)=rk;
sumrkA(:,i)=sum(rk>0);
end
sumrkA

