function [x, iter] = FqpEicp(A, B, x0, maxIt, maxItflqp, maxItsub, strategy, eps, epsflqp, epsbbp)
x = x0;
[~, n] = size(A);
iter = 0;
while 1   
x = FLQP(x0, A, B, n, strategy, maxItflqp, maxItsub, epsflqp, epsbbp);
e = norm(x - x0);
if e < eps || iter > maxIt
    break;
else
    x0 = x;
%     lambda = (x' * A * x) / (x' * B * x)
end
iter = iter + 1;
end
