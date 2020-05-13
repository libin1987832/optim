function [ xk,res] = splitD(Q,d,iter);
[m,n]=size(Q)
Di=diag(Q);
Di(Di<0)=0;
Di=diag(1./Di);
tol=1e-12;
x0=zeros(n,1);
for i=1:iter
    xk=x0-Di*(d+Q*x0);
    xk(xk<0)=0;
    x0=xk;
    res=test_valid(Q,d,xk);
    if res<tol
        break;
    end
end