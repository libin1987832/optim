clear
clc

m = 30000;
n = 2000;

A = 2 * rand(m , n)-1;

b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
t=clock;
xk=krylovk(A,b,80);
tf_ky=etime(clock,t);
res_ky = norm(A'*(A*xk-b));
fprintf('krylov residual %g,time %g.\n',res_ky,tf_ky);
maxit = 50000;
tol = [];
exactx = [];
debug = 0;
t=clock;
[x,iter,error_k,iter_k,index_k] = randomizedGaussSeidel(A, b, x0,maxit,tol,exactx,debug);
tf_gs=etime(clock,t);
res_gs=norm(A'*(A*x-b));
fprintf('rand residual %g,time %g.\n',res_gs,tf_gs);
p=2;
t=clock;
[x,iter,error_k,iter_k,index_k] = wrandomizedGaussSeidel(A, b, x0,p,maxit,tol,exactx,debug);
tf_wgso=etime(clock,t);
res_wgso=norm(A'*(A*x-b));
fprintf('weight rand residual %g,time %g.\n',res_wgso,tf_wgso);
% t=clock;
% opts.strategy=3;
% opts.p=2;
% opts.xstar = xs;
% [x,Out]=dARGauss_Seidel(A,b,opts);
% tf_wgs=etime(clock,t);
% res_wgs=norm(A*x-b);
% fprintf('dAR rand residual %g,time %g.\n',res_wgs,tf_wgs);

% figure
% display=1:1:1000;
% h=semilogy(display, Out.error(display), 'k.');