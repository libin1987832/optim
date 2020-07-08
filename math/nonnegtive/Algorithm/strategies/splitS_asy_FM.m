function [xkA,rpk,cmax] = splitS_asy_FM(A,b,Q,R,x0,iter,cmax,r0)
[m,n]=size(A);
if iter<1
    xk=x0;
    xkA=xk;
else
    xkA=zeros(n,iter);
%    tol=1e-12;
    % x0=zeros(n,1);
    for i=1:iter
        [xk,rpk]=FixedM(x0,Q,R,A,b,r0);
        r0=rpk;
        if i>1
            d1=norm(xk-x0);
            d2=norm(x0-xn1);
            c=d1/d2;
            if c>cmax
                cmax=c;
            end
        end
        xn1=x0;
        x0=xk;
        xkA(:,i)=xk;
    end
end
