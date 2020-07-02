function xk=MPRGP(A,b,x0,L,a,delta,maxIter)
k=0;r=A*x0-b;
tol=1e-15;
[m,n]=size(r);
[fx,bx]=fbxf(x0,m,tol,r);
p=fx;
bxn=bx'*bx;
while bxn+fx'*fx>delta
    f1x=(x0-l)./a;
    fxc=find(f1x>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx
       sqrt()
       disp(['conject gradient:',num2str(r'*r)]);
    else
       FM() 
       disp(['go to implement space:',num2str(r'*r)]);
    end
    k=k+1;
    bxn=bx'*bx;
    x0=xk;
end
end

function [amax,xc,res]=computeMaxStep(x0,d)
grt0=(d>0);
% if all indredient great 0 then natural sartisfied
if sum(grt0)<1
    amax=100;
    res=0;
else
    amax=min(x0(grt0)./d(grt0));
    res=1;
end
xc=x0+amax*d;
end

function [f,b]=fbxf(x0,m,tol,p)
AA=find(x0<tol);
FF=setdiff(1:m,AA);
b=zeros(m,1);
f=zeros(m,1);
b(AA)=p(AA);
f(FF)=p(FF);
b(b>0)=0;
end


function f=fxf(x0,m,tol,p)
FF=find(x0>tol);
f=zeros(m,1);
f(FF)=p(FF);
end

function b=bxf(x0,m,tol,p)
AA=find(x0<tol);
b=zeros(m,1);
b(AA)=p(AA);
b(b>0)=0;
end
