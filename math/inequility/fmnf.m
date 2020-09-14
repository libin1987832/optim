function [xk, rpk]=fmnf(A,b,xk,Q,rpk,nf)
for i = 1:nf
    r = rpk;
    r(r<0) = 0;
    Qr = Q' * r;
    % compute min increase
    uk = R \ Qr;
    xk = xk+uk;
    rpk = b - A * xk;
end