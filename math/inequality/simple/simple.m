% parameter: A b x0 n obvious
% QR is the decomposition for A if Dax method
% R may be steplength if the gradient method
% rpk is b-Ax
% compuataion: GHA: A' * r;,R * d,A * xk FM Q' * r  R \ Qr A * xk
function [xkA,rpk] = simple(A,b,xk,n,rpk,nf,D)
xkA = zeros(n,nf);
for i = 1 : nf
%    [xk,rpk] = FMGS2(xk,A,b,D,rpk,1);
%    rpk = b - A * xk;
     rpk(rpk<0)=0;
<<<<<<< HEAD
      uk=krylovk(A,rpk,10);
=======
      uk=krylovk(A,rpk,50);
>>>>>>> 68f37df8bd2b2fc3e8ce82a94bceedf6c5f0fdf4
      xk = xk+uk;
      rpk = b-A*xk;
    xkA(:,i) = xk;
end
