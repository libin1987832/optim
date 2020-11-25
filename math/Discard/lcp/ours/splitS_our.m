function xkA= splitS_our(Q,d,s,x0,iter)
xkA=[];
if iter<1
    xk=x0;
    xkA=[xkA xk];
else
    addpath('../test')
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
        xn1=x0;
        x0=xk;
        xkA=[xkA xk];
    end
end
