function [res, xmin, Armax] = kktResidual(A, b, x)
r = b - A * x;
r(r<0) = 0;
Ar = -A' * r;
mxAr = min(x,Ar);
res = sqrt(mxAr' * mxAr);
xmin = min(x);
Armax = max(abs(Ar));