function [x, error] = sqpEicp(A, B, M, x0, eps, maxIT, debug)
error = 0;
x = x0;
if debug
    error = zeros(1, maxIT);
end
k=0;
while 1
d = searchdir(A, B, M, x);
if norm( d ) < eps || k > maxIT
    break;
end
alpha = searchstep(A, B, M, d);
x = x + alpha * d;
k = k + 1;
if debug
    error(1,k) = 0.5 * x' * A * x;
end
end