function [rpk, normr, xmin, Ar, normKKT , face1, face2] = kktResidual(A, b, x , rpk, type)
% if the parametre is greater than 3 , save the computation of  rpk
if nargin < 4 || isempty(rpk)
rpk = b - A * x;
end
r = rpk;
r(r<0) = 0;
normr = sqrt( r' * r );
%xmin = min( x );
xmin = norm( x );
Ar = -1;
normKKT = -1;
% if the parametre is greater than 4 , save the computation of  Ar 
if nargin == 5
    Ar = -A' * r;
   
    mxAr = min( x , Ar );
    normKKT = sqrt( mxAr' * mxAr );
    %normKKT = max(abs(x.*Ar));
    Ar = min(Ar);
    face1 = sum(rpk>-1e-15);
    face2 = sum(x>1e-10);
end