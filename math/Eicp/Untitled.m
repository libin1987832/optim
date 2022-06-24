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
d =  Ax0 / a;
Ax0 = A * x0;
Bx0 = B * x0;
xBx = 0.5 * x0' * Bx0;
xAx = 0.5 * x0' * Ax0;
a = xAx / xBx;
eps = 1e-5;
ninf1 = 2 * n;
maxIt = 1000;
[x1,fval,exitflag,output,lambda] = quadprog(B, -Ax0 / a, -Ax0', -xAx , ones(1, n), 1, zeros(n, 1), [ ]);
[x0, F, iter, ninf, testwx] = BBP2(Ax0, B, a, xAx, maxIt, ninf1, eps, 1, 0);

[x, F, iter, ninf, testwx] = BBP2(Ax0, B, a, xAx, maxIt, ninf1, eps, 0, 0);
M = [ B -ones( n , 1 ) -Ax0; ones( 1 , n ) 0 0; Ax0' 0 0 ]; 
M * [x1;-lambda.eqlin;-lambda.ineqlin] + [-Ax0/a;-1;-xAx]
[x0,x,x1]