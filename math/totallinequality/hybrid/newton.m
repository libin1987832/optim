function [xs,rpk] = newton(A,b,n,rpk,x0)
tol = 1e-15;
AA = (rpk > -tol);
RR = (x0 > tol);
% subspace
AI = A(AA, RR);
bI = rpk( AA );
u = lsqminnorm(AI, bI);
p = zeros(n, 1);
p(RR) = u;
% debug for test.m nf = 2
[alpha, aranges, retcode] = arraySpiece(A, b, x0, p);
xs = x0 + alpha * p;
rpk = b - A * xs;
