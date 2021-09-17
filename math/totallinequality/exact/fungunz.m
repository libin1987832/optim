function [f,g]=fungunz(x,A,b,AtA,Ab)
f = norm( A*x - b)^2;
g = 2*( AtA*x - Ab );