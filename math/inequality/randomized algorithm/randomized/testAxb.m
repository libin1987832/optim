clear
clc

m = 10000;
n = 1000;

A = 2 * randn(m , n)-1;
AA=A.*A;
columua=sum(AA,1);
A = AA./repmat(columua,m,1);
b = 2 * randn(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
xs = A\b;
norm(A*xs-b)
xk=krylovk(A,b,2);
norm(A*xk-b)
opts.strategy=2;
opts.p=2;
opts.xstar = xs;
[x,Out]=dARGauss_Seidel(A,b,opts);
norm(A*x-b)
figure
display=1000:100:Out.iter;
h=semilogy(display, Out.error(display), 'k.');