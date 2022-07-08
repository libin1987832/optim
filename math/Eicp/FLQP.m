function [x, iter, Abpp,error] = FLQP(a, B, Ax0, xAx, n, strategy, maxIt, maxItsub, eps, epsbbp, debug)
%Ax0 = A * x0;
%Bx0 = B * x0;
%xBx = x0' * Bx0;
%xAx = x0' * Ax0;
%a = (xAx) / xBx;
%a = xBx / xAx;
a0 = a + 10;
iter = 0;
ninf1 = 2 * n;
%fun = @(x) -(2 * Ax0' * x - xAx) / (x' * B * x);
%opts = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
%x = fmincon(fun,x0,-Ax0', -xAx ,ones(1, n), 1, zeros(n, 1), inf(n,1));
if debug
fun = @(x) (x'* Ax0)/(x' * B * x);
error=zeros(1,maxIt);
end
F = 1:n;
nitBBs = 0;
 while abs(a - a0) > eps && abs(a) > eps && iter < maxIt
a0 = a;
iter = iter + 1;
if strategy == 1
    [x, F, ~, ninf, testwx] = BBP2(Ax0, B, n, a, 0.5 * xAx, maxItsub, ninf1, epsbbp, 0, debug);
    Bx0 = B * x;
    xBx = x' * Bx0;
    a = xBx / ( 2 * Ax0' * x - xAx);
elseif strategy == 2
     [x, ~, F, nitBB, ~] = BBP3(a * B, -Ax0, F, maxItsub, epsbbp, 0);
     a = (Ax0' * x) / (x' * B * x);
     nitBBs = nitBBs + nitBB;
else
    opts = optimset('Display','off');
    [x,fval,exitflag,output,lambda] = quadprog(B, -a * Ax0, -Ax0', -0.5*xAx , ones(1, n), 1, zeros(n, 1), [ ], [], opts);
    Bx0 = B * x;
    xBx = x' * Bx0;
    a = xBx / ( 2 * Ax0' * x - xAx);
 %   x = fmincon(fun,x0,-Ax0', -xAx ,ones(1, n), 1, zeros(n, 1), inf(n,1));
end
if debug
error(iter) = fun(x);
end

 %a  = ( 2 * Ax0' * x0 - xAx) / xBx;
 end
Abpp = nitBBs / iter;
% Ax0 = A * x0;
% xAx = 0.5 *x0' * Ax0;
% d =  Ax0 / a;
% eps = 1e-5;
% ninf1 = 2 * n;
% maxIt = 1000;

