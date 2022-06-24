function x = FLQP(x0, A, B, n, strategy, eps)
Ax0 = A * x0;
Bx0 = B * x0;
xBx = x0' * Bx0;
xAx = x0' * Ax0;
a = - (xAx - Ax0' * x0) / xBx;
a0 = a + 10;
maxIt = 1000;
epssub = 1e-5;
ninf1 = 2 * n;
fun = @(x) -(2 * Ax0' * x - xAx) / (x' * B * x);
%opts = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
x = fmincon(fun,x0,-Ax0', -xAx ,ones(1, n), 1, zeros(n, 1), inf(n,1));
% while abs(a - a0) > eps && abs(a) > eps
% a0 = a;
% if strategy == 1
%     [x, F, iter, ninf, testwx] = BBP2(Ax0, B, n, a, xAx, maxIt, ninf1, epssub, 1, 0);
% else
%   %  [x,fval,exitflag,output,lambda] = quadprog(a * B, -Ax0, -Ax0', -xAx , ones(1, n), 1, zeros(n, 1), [ ]);
%     x = fmincon(fun,x0,-Ax0', -xAx ,ones(1, n), 1, zeros(n, 1), inf(n,1));
% end
% % Ax0 = A * x;
%  Bx0 = B * x;
%  xBx =  * x' * Bx0;
% % xAx = 0.5 * x' * Ax0;
% a = - (xAx - Ax0' * x) / xBx;
% end
% Ax0 = A * x0;
% xAx = 0.5 *x0' * Ax0;
% d =  Ax0 / a;
% eps = 1e-5;
% ninf1 = 2 * n;
% maxIt = 1000;

