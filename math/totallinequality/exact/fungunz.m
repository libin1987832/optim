function [f,g]=fungun(x,A,b)
f = norm( A*x - b)^2;
g = 2*A'*(A*x-b);