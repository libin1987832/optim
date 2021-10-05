% parameter: A b x0 n obvious
% QR is the decomposition for A if Dax method
% R may be steplength if the gradient method
% rpk is b-Ax
% compuataion: GHA: A' * r;,R * d,A * xk FM Q' * r  R \ Qr A * xk
function [xkA,rpk,resA] = simple(A,b,xk,n,Q,R,rpk,nf,type)
xkA = zeros(n,nf);
resA = zeros(1,nf);
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
           % uk = R \ Qr;
           uk=zeros(n,1);
           for k=n:-1:1
               for j=k:n
                   Qr(k)=Qr(k)-R(k,j)*uk(j);
               end
               uk(k)=Qr(k)/R(k,k);
           end
           xk = xk+uk;
           rpk = b - A * xk;
           xkA(:,i) = xk;
           [~, r, normr, normAr] = residual(A,b,xk,rpk);
           resA(1, i) = normAr;
        end        
end