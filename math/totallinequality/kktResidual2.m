% ||Ax-b|| 
function [rpk, normr, xmin, Ar, normKKT , faceX, faceA] = kktResidual2(A, b, x , rpk, type)
% if the parametre is greater than 3 , save the computation of  rpk
if nargin < 4 || isempty(rpk)
rpk = A * x - b;
end
r = rpk;
r(r<0) = 0;
normr = 0.5 * norm(r,2)^2;
%xmin = min( x );
xmin = norm( x );
Ar = -1;
normKKT = -1;
% if the parametre is greater than 4 , save the computation of  Ar 
if nargin == 5
    Ar = A' * r;
    mxAr = min( x , Ar );
    %normKKT = sqrt( mxAr' * mxAr );
    normKKT = max(abs(x.*Ar));
    %Ar = min(Ar);
    faceA = sum(rpk>-1e-15);
    faceX = sum(x>1e-10);
end