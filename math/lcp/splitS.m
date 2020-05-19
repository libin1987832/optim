function [xk,res] = splitS(Q,d,s,x0,iter)
[m,n]=size(Q);
Di=diag(Q);
Di(Di<0)=0;
% Di=1./Di;
tol=1e-12;
% x0=zeros(n,1);

for i=1:iter
    xk=zeros(n,1);

    for j=1:n
%         sd=s*Di(j)*(d(j)+Q(j,:)*[xk(1:(j-1));x0(j:n)]);
        sd=s*(d(j)+Q(j,:)*[xk(1:(j-1));x0(j:n)])/Di(j);
        if x0(j)>sd
            xk(j)=x0(j)-sd;
        end
    end
%     xk(xk<0)=0;
    x0=xk;
    res=test_valid(Q,d,xk);
    if res<tol
        break;
    end
end