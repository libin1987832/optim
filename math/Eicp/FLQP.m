function [x, iter] = FLQP(x0, A, B, n, strategy, maxIt, maxItsub, eps, epsbbp)
Ax0 = A * x0;
Bx0 = B * x0;
xBx = x0' * Bx0;
xAx = x0' * Ax0;
%a = (xAx) / xBx;
a = xBx / xAx;
a0 = a + 10;
iter = 0;
ninf1 = 2 * n;
%fun = @(x) -(2 * Ax0' * x - xAx) / (x' * B * x);
%opts = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
%x = fmincon(fun,x0,-Ax0', -xAx ,ones(1, n), 1, zeros(n, 1), inf(n,1));
 while abs(a - a0) > eps && abs(a) > eps && iter < maxIt
a0 = a;
iter = iter + 1;
if strategy == 1
    [x, F, ~, ninf, testwx] = BBP2(Ax0, B, n, a, 0.5 * xAx, maxItsub, ninf1, epsbbp, 0, 0);
else
    opts = optimset('Display','off');
    [x,fval,exitflag,output,lambda] = quadprog(B, -a * Ax0, -Ax0', -0.5*xAx , ones(1, n), 1, zeros(n, 1), [ ], [], opts);
 %   x = fmincon(fun,x0,-Ax0', -xAx ,ones(1, n), 1, zeros(n, 1), inf(n,1));
end
 Bx0 = B * x;
 xBx = x' * Bx0;
 a = xBx / ( 2 * Ax0' * x - xAx);
 %a  = ( 2 * Ax0' * x0 - xAx) / xBx;
end
% Ax0 = A * x0;
% xAx = 0.5 *x0' * Ax0;
% d =  Ax0 / a;
% eps = 1e-5;
% ninf1 = 2 * n;
% maxIt = 1000;

