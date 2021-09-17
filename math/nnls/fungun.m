function [f,g]=fungun(A,b,x)
f = norm( A*x - b)^2;
g = 2*A'*(A*x-b);