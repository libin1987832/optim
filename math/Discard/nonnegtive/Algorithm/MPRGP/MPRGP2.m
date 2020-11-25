function xk=MPRGP2(A,b,x0,L,a,delta,maxIter)
k=0;r=A*x0-b;
z0=r;
z0(z0<0)=0;
tol=1e-15;
[m,n]=size(r);
[fx,bx]=fbxf(x0,m,tol,r);
p=fx;
bxn=bx'*bx;

AA=find(z0<tol);
FF=setdiff(1:m,AA);
AI=A(AA,:);
rI=AI*x0-bI;
p=rI;
while bxn+fx'*fx>delta
    f1x=x0./a;
    fxc=find(f1x>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx 
        acg=(rI'*p)/(p'*AI*p);y=x0-acg*p;
        azf=z(FF)./(AI*p);
        
        af=computeMaxStep(x0,p);
        if acg<af
            xk=y;r=r-acg*A*p;[fy,bx]=fbxf(y,m,l,tol,r);
            grma=(fy'*A*p)/(p'*A*p);
            p=fy-grma*p;
            fx=fy;
            disp(['conject gradient:',num2str(r'*r)]);
        else
            xk2=x0-af*p;r=r-af*A*p;
            fx2=fxf(xk2,m,l,tol,r);
            xk=xk2-a*fx2;
            xk(xk<l)=l;
            r=A*xk-b;[p,bx]=fbxf(xk,m,l,tol,r);
            fx=p;
            disp(['out range in the subspace:',num2str(r'*r)]);
        end
    else
        d=bxf(x0,m,l,tol,r);
        acg=(r'*d)/(d'*A*d);
        xk=x0-acg*d;r=r-acg*A*d;[p,bx]=fbxf(xk,m,l,tol,r);
        fx=p;
        disp(['go to implement space:',num2str(r'*r)]);
    end
    k=k+1;
    bxn=bx'*bx;
    x0=xk;
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

function [f,b]=fbxf(x0,m,l,tol,p)
AA=find(x0<l+tol);
FF=setdiff(1:m,AA);
b=zeros(m,1);
f=zeros(m,1);
b(AA)=p(AA);
f(FF)=p(FF);
b(b>0)=0;
end


function f=fxf(x0,m,l,tol,p)
FF=find(x0>l+tol);
f=zeros(m,1);
f(FF)=p(FF);
end

function b=bxf(x0,m,l,tol,p)
AA=find(x0<l+tol);
b=zeros(m,1);
b(AA)=p(AA);
b(b>0)=0;
end
