function [x,rpk] = projectfixed(A,b,x0,rpk,alpha)

rpk(rpk<0) = 0;
g = A' * rpk;
x = x0 + alpha * g;
x(x<0) = 0;

rpk = b - A * x;