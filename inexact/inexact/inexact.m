function [x0,f1]=inexact(b,A,x0,belt,u)
[msize,nsize]=size(x0);
while 1
    y0=b-A*x0;
    y0(y0<=0)=0;
    
    b0=y0+A*x0;
    s0=A'*y0;
    
    sN=s0;
    alph=sN'*sN/((A*sN)'*(A*sN));
    mk=0;
    lamd=diag(ones(msize,1));
%     s0=-1*s0;
    for i=1:msize
        if 0<x0(i) & s0(i)<0
            lamd(i,i)=-1*x0(i)/s0(i);
        end
    end
    while 1
        x1=x0+belt^mk*lamd*s0;
        x1(x1<0)=0;
        left=(b0-A*x1)'*(b0-A*x1);
        right=y0'*y0-2*u*s0'*(x1-x0);
        if left<=right
            break;
        end
        mk=mk+1;
    end
    f0=b-A*x0;
    f0(f0<0)=0;
    f0=0.5*(f0'*f0);
    f1=b-A*x1;
    f1(f1<0)=0;
    f1=0.5*(f1'*f1);
    fprintf('The mk is %f; x0,x1:%f,%f; f0,f1:%f,%f\n',mk,x0,x1,f0,f1);
    x0=x1;
    if (f1-f0)'*(f1-f0)<0.00001
        break;
    end
end