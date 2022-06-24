function x = FqpEicp(A, B, x0, eps)
x = x0;
[~, n] = size(A);
while 1   
x = FLQP(x0, A, B, n, 0, 1e-10);
e = norm(x - x0);
if e < eps
    break;
else
    x0 = x;
end
end
