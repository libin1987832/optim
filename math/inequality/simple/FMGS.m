function [xk,rk] = FMGS(x0,A,b,D,rk,niter)
[m,n]=size(A);
xk=x0;
zk=-rk;
zk(zk<0)=0;
% compute min increase
for j=1:niter
    for i=1:n
        r=rk+zk;
        xk(i) = xk(i) + (A(:,i)'*r)/D(i);
        rk = b - A*xk;
    end
end