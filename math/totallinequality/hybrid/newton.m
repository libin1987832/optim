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
xa = [0.019:0.00001:0.021];
ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
 pxy={};
% pxy(1).X = knot;
% pxy(1).Y = knoty;
pxy(1).X = xa;
pxy(1).Y = ya;
figure
p1 = arrayfun(@(a) plot(a.X,a.Y),pxy);
hold on