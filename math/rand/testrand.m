m = 2000;
n = 500;
A = rand(m, n);
b = rand(m, 1);
x0 = zeros(n, 1);
tol = 1e-13;
maxit = 1000;
[x0, resvec, randp] = randkaczmarz(A, b, x0, tol, maxit);
