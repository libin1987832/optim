% parameter: A b x0 n obvious
% QR is the decomposition for A if Dax method
% R may be steplength if the gradient method
% rpk is b-Ax
% compuataion: GHA: A' * r;,R * d,A * xk FM Q' * r  R \ Qr A * xk
function [xkA,rpk] = simple(A,b,xk,n,rpk,nf,D)
xkA = zeros(n,nf);
for i = 1 : nf
    xk = FMGS(xk,A,b,D,rpk,1);
    rpk = b - A * xk;
    xkA(:,i) = xk;
end