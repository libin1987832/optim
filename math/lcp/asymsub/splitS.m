function [xk,cmax] = splitS(Q,d,s,x0,iter,cmax)
addpath('./test')
[m,n]=size(Q);
Di=diag(Q);
Di(Di<0)=0;
tol=1e-12;
% x0=zeros(n,1);
for i=1:iter
    xk=zeros(n,1);
    for j=1:n
        sd=s*(d(j)+Q(j,:)*[xk(1:(j-1));x0(j:n)])/Di(j);
        if x0(j)>sd
            xk(j)=x0(j)-sd;
        end
    end
    if i>2
        d1=norm(xk-x0);
        d2=norm(x0-xn1);
        c=d1/d2;
        if c>cmax
           cmax=c;
        end
    end
    xn1=x0;
    x0=xk;
end