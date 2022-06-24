function x = FLQP(x0, A, B, a, n, strategy)
Ax0 = A * x0;
xAx = x0' * Ax0;
d =  Ax0 / a;
if strategy == 1
    x = BBP2(A, B, b, F, x0, maxIt, ninf1, eps, strategy, debug);
else
    [x, fval, exitflag] = quadprog(-B, d, -d', - 0.5 * xAx , ones(1, n), 1, zeros(n, 1), [ ]);
end