function [x, iter, italg1] = FqpEicp(A, B, x0, maxIt, maxItflqp, maxItsub, strategy, eps, epsflqp, epsbbp)
x = x0;
[~, n] = size(A);
iter = 0;
itflqps = 0;
while 1  
iter = iter + 1;
[x, itflqp] = FLQP(x0, A, B, n, strategy, maxItflqp, maxItsub, epsflqp, epsbbp);
itflqps = itflqps + itflqp;
e = norm(x - x0);
if e < eps || iter > maxIt-1
    break;
else
    x0 = x;
%     lambda = (x' * A * x) / (x' * B * x)
end
end
italg1 = itflqps / iter;
