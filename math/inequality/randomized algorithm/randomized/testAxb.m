clear
clc

m = 1000;
n = 1000;

A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
xs = A\b;
xk=krylovk(A,b,10);
norm(A*xk-b)
opts.strategy=2;
opts.p=20;
opts.xstar = xs;
[x,Out]=dARGauss_Seidel(A,b,opts);
norm(A*x-b)
figure
h=semilogy(Out.iter, Out.error, 'k.');