% parameter: A b x0 n obvious
% QR is the decomposition for A if Dax method
% R may be steplength if the gradient method
% rpk is b-Ax
% compuataion: GHA: A' * r;,R * d,A * xk FM Q' * r  R \ Qr A * xk
% function [xkA,rpk] = simpleNqr(A,b,xk,n,steplengthOrk,rpk,nf,type)
function [xkA,rpk] = simpleNqr(A,b,xk,n,rpk,param)
nf = param.nf;
type = param.type;
steplengthOrk = param.steplengthOrk;
xkA = zeros(n,nf);
switch upper(type)
    case 'GHA'
        for i = 1 : nf
            r = rpk;
            r(r<0) = 0;
            d = A' * r;
            xk = xk + steplength * d;
            rpk = b - A * xk;
            xkA(:,i) = xk;
        end                
    otherwise 
        for i = 1 : nf
            r = rpk;
            r(r<0) = 0;
            uk=krylovk(A,r,steplengthOrk);
            xk = xk+uk;
            rpk = b - A * xk;
            xkA(:,i) = xk;
        end        
end