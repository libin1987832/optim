function x = FqpEicp(A, B, x0, eps)
x = x0;
while 1   
x = FLQP(x0, A, B, n, 1);
e = norm(x - x0);
if e < eps
    break;
else
    x0 = x;
end
end
