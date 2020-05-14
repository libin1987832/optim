function [xk,res] = splitS(Q,d,s,x0,iter)
[m,n]=size(Q)
Di=diag(Q);
Di(Di<0)=0;
Di=1./Di;
tol=1e-12;
% x0=zeros(n,1);
xk=x0;
for i=1:iter
    xk(1)=x0(1)-s*Di(1)*(d(1)+Q(1,:)*x0);
    for j=2:n
        sd=s*Di(j)*(d(j)+Q(j,1:(j-1))*xk(1:(j-1))+Q(j,j:n)*x0(j:n));
        xk(j)=x0(j)-sd;
    end
    xk(xk<0)=0;
    x0=xk;
    res=test_valid(Q,d,xk);
    if res<tol
        break;
    end

end