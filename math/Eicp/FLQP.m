function x = FLQP(x0, A, B, n, strategy)
Ax0 = A * x0;
Bx0 = B * x0;
xBx = 0.5 * x0' * Bx0;
xAx = 0.5 * x0' * Ax0;
a = xAx / xBx;
eps = 1e-5;
ninf1 = 2 * n;
maxIt = 1000;
if strategy == 1
    x = BBP2(Ax0, B, a, xAx, maxIt, ninf1, eps, 1, 0);
else
    [x,fval,exitflag,output,lambda] = quadprog(B, -Ax / a, -Ax', -xAx , ones(1, n), 1, zeros(n, 1), [ ]);
end
% Ax0 = A * x0;
% xAx = 0.5 *x0' * Ax0;
% d =  Ax0 / a;
% eps = 1e-5;
% ninf1 = 2 * n;
% maxIt = 1000;

