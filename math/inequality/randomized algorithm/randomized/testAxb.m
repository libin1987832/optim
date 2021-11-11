clear
clc

m = 10000;
n = 1000;

A = 2 * randn(m , n)-1;
% AA=A.*A;
% columua=sum(AA,1);
% A = AA./repmat(columua,m,1);
b = 2 * randn(m , 1)-1;
% b=A*ones(n,1);
x0 = zeros(n , 1);
xs = A\b;
norm(A*xs-b)
xk=krylovk(A,b,2);
norm(A*xk-b)
maxit = 1000;
tol = [];
exactx = [];
debug = 0;
[x,iter,error_k,iter_k,index_k] = randomizedGaussSeidel(A, b, x0,maxit,tol,exactx,debug);
norm(A*x-b)
p=2;
[x,iter,error_k,iter_k,index_k] = wrandomizedGaussSeidel(A, b, x0,p,maxit,tol,exactx,debug);
norm(A*x-b)
opts.strategy=2;
opts.p=2;
opts.xstar = xs;
[x,Out]=dARGauss_Seidel(A,b,opts);
norm(A*x-b)
figure
display=1:10:1000;
h=semilogy(display, Out.error(display), 'k.');