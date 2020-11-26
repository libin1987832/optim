function [r, Ar, xmin, normKKT] = kktResidual(A, b, x)
r = b - A * x;
r(r<0) = 0;
Ar = -A' * r;
mxAr = min(x,Ar);
normKKT = sqrt(mxAr' * mxAr);
xmin = min(x);
