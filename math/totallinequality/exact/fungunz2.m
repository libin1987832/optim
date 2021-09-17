function [f,g]=fungunz2(x,A,b,AtA,Ab)
f = norm( A*x - b)^2;
g =  2*A'*(A*x-b);