function xk=MPRGP(A,b,x0,L,a,delta,l)
k=0;r=A*x0-b;
tol=1e-15;
[m,n]=size(r);
fx=fxf(x0,m,l,tol,r);
bx=bxf(x0,m,l,tol,r);
p=fx;
bxn=bx'*bx;
fxn=fx'*fx;
while bxn+fxn>delta
    f1x=(x0-l)./a;
    fxc=find(flx>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx
        acg=(r'*p)/(p'*A*p);y=x0-acg*p;
        af=computeMaxStep(x0,p,l);
        if acg<af
            xk=y;r=r-acg*A*p;fy=fxf(y,m,l,tol,r);
            grma=(fy'*A*p)/(p'*A*p);
            p=fy-grma*p;
            fx=fy;
        else
            xk2=x0-af*p;r=r-af*A*p;
            fx2=fxf(xk2,m,l,tol,r);
            xk=xk2-a*fx2;
            xk(xk<l)=l;
            r=A*xk-b;p=fxf(xk,m,l,tol,r);
            fx=p;
        end
    else
        d=bxf(x0,m,N,l,tol,r);
        acg=(r'*d)/(d'*A*d);
        xk=x0-acg*d;r=r-acg*A*d;p=fxf(xk,m,l,tol,r);
        fx=p;
    end
    k=k+1;
    bxn=bx'*bx;
    fxn=fx'*fx;
end
end

function [amax,xc,res]=computeMaxStep(x0,d,l)
grt0=(d>0);
% if all indredient great 0 then natural sartisfied
if sum(grt0)<1
    amax=100;
    res=0;
else
    amax=min((x0(grt0)-l)./d(grt0));
    res=1;
end
xc=x0+amax*d;
end

function f=fxf(x0,m,l,tol,p)
FF=find(x0>l+tol);
f=zeors(m,1);
f(FF)=p(FF);
end

function b=bxf(x0,m,l,tol,p)
AA=find(x0<l+tol);
b=zeors(m,1);
b(AA)=p(AA);
b(b>0)=0;
end

function [f,b]=fbxf(x0,m,l,tol,p)
AA=find(x0<l+tol);
FF=setdiff(1:N,AA);
b=zeors(m,1);
f=zeors(m,1);
b(AA)=p(AA);
f(FF)=p(FF);
b(b>0)=0;
end