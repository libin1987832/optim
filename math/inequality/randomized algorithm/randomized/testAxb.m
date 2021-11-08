clear
clc

m = 1000;
n = 1000;

A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);

 xk=krylovk(A,b,10);
norm(A*xk-b)
 [x,Out]=dARGauss_Seidel(A,b,[]);
 norm(A*x-b)