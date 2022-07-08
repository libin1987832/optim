function [x, iter, italg1] = FqpEicp(A, B, x0, maxIt, maxItflqp, maxItsub, strategy, eps, epsflqp, epsbbp)
x = x0;
[~, n] = size(A);
iter = 0;
itflqps = 0;
Ax0 = A * x0;
Bx0 = B * x0;
xBx = x0' * Bx0;
xAx = x0' * Ax0;
%a = (xAx) / xBx;
a = xBx / xAx;
while 1  
iter = iter + 1;
%[x, itflqp] = FLQP(x0, A, B, n, strategy, maxItflqp, maxItsub, epsflqp, epsbbp);
[x, itflqp] = FLQP(a, B, Ax0, xAx, n, strategy, maxIt, maxItsub, eps, epsbbp);
itflqps = itflqps + itflqp;
e = norm(x - x0);
if e < eps || iter > maxIt-1
    break;
end
x0 = x;
Ax0 = A * x0;
Bx0 = B * x0;
xBx = x0' * Bx0;
xAx = x0' * Ax0;
a = xBx / xAx;
ai = 1./a;
BAx = ai * Bx0 - Ax0;
if e < eps || iter > maxIt-1
    break;
end
end
italg1 = itflqps / iter;
