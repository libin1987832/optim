% x0 initial value  b
% xk out value 
%r0=(b-Ax0) rk=(b-Axk)  
function [xk,rk] = FixedGS(x0,A,b,D,niter)
[m,n]=size(A);
xk=x0;
% compute min increase
for j=1:niter
    Axk = A*xk; 
    zk=Axk-b;
    zk(zk<0)=0;
    bz=b+zk;
    for i=1:n
        r=bz-Axk;
        xk(i)=xk(i)-(A(:,i)'*r)/D(i);
        Axk = A*xk;
    end
end
rk = b-Axk;