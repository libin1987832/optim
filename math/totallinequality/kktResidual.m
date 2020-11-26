function [r, normr, xmin, Ar, normKKT] = kktResidual(A, b, x , type)
r = b - A * x;
r(r<0) = 0;
normr = sqrt( r' * r );
xmin = min( x );
Ar = -1;
normKKT = -1;
if nargin == 4
    Ar = -A' * r;
    mxAr = min( x , Ar );
    normKKT = sqrt( mxAr' * mxAr );
end
