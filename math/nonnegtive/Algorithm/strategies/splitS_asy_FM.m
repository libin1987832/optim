function [xkA,rk,cmax] = splitS_asy_FM(A,b,Q,R,x0,iter,cmax)
xkA=[];
if iter<1
    xk=x0;
    xkA=[xkA xk];
else
    [m,n]=size(A);
%    tol=1e-12;
    % x0=zeros(n,1);
    for i=1:iter
        [xk,r0,rk]=FixedM(x0,Q,R,A,b);
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
        xkA=[xkA xk];
    end
end
