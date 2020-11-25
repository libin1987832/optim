function [xkA, rpk]=fmnf(A,b,xk,n,Q,R,rpk,nf)
xkA = zeros(n,nf);
for i = 1:nf
    r = rpk;
    r(r<0) = 0;
    Qr = Q' * r;
    % compute min increase
    uk = R \ Qr;
    xk = xk+uk;
    rpk = b - A * xk;
    xkA(:,i) = xk;
end