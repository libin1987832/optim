function [xkA,rpkA] = simpleQR(A,b,xk,n,m,rpk,nf,Q,R)
xkA = zeros(n,nf);
rpkA = zeros(m,nf);
for i = 1 : nf
   [xk,rpk]=FMQR(xk,Q,R,A,b,rpk);
    xkA(:,i) = xk;
    rpkA(:,i) = rpk;
end