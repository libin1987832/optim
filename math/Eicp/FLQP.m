function x = FLQP(x0, A, B, a, n, strategy)
Ax0 = A * x0;
xAx = 0.5 *x0' * Ax0;
d =  Ax0 / a;
eps = 1e-5;
ninf1 = 2 * n;
if strategy == 1
    x = BBP2(d, B, xAx, maxIt, ninf1, eps, 0, 0);
else
    [x, fval, exitflag] = quadprog(-B, d, -d', -xAx , ones(1, n), 1, zeros(n, 1), [ ]);
end