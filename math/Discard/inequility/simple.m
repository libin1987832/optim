% parameter: A b x0 n obvious
% QR is the decomposition for A if Dax method
% R may be steplength if the gradient method
% rpk is b-Ax
function [xkA,rpk] = simple(A,b,xk,n,Q,R,rpk,nf,type)
xkA = zeros(n,nf);
switch upper(type)
    case 'GHA'
        for i = 1 : nf
            r = rpk;
            r(r<0) = 0;
            d = A' * r;
            xk = xk + R * d;
            rpk = b - A * xk;
            xkA(:,i) = xk;
        end                
    otherwise 
        for i = 1 : nf
            r = rpk;
            r(r<0) = 0;
            Qr = Q' * r;
            % compute min increase
            uk = R \ Qr;
            xk = xk+uk;
            rpk = b - A * xk;
            xkA(:,i) = xk;
        end        
end