function [x, iter] = FqpEicp(A, B, x0, maxIt, strategy, eps)
x = x0;
[~, n] = size(A);
iter = 0;
while 1   
x = FLQP(x0, A, B, n, strategy, 1e-15);
e = norm(x - x0);
if e < eps || iter > maxIt
    break;
else
    x0 = x;
%     lambda = (x' * A * x) / (x' * B * x)
end
iter = iter + 1;
end
