N = 5;
n = N;
D = diag(rand(N,1));
U = orth(rand(N,N));
A = U' * D * U;
D = diag(rand(N,1));
U = orth(rand(N,N));
B = U' * D * U;
x0 = rand(N, 1);
x0 = x0 ./ sum(x0);
a = 0.5;
Ax0 = A * x0;
xAx = 0.5 *x0' * Ax0;
d =  Ax0 / a;
Ax0 = A * x0;
Bx0 = B * x0;
xBx = 0.5 * x0' * Bx0;
xAx = 0.5 * x0' * Ax0;
a = xAx / xBx;
eps = 1e-5;
ninf1 = 2 * n;
maxIt = 1000;
[x, F, iter, ninf, testwx] = BBP2(Ax0, B, a, xAx, maxIt, ninf1, eps, 1, 0);
[x1,fval,exitflag,output,lambda] = quadprog(B, -Ax0 / a, -Ax0', -xAx , ones(1, n), 1, zeros(n, 1), [ ]);

[x,x1]