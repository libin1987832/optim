function MPRGP(A,b,x0,L,a,delta,l)
k=0;r=A*x0-b;p=r;
tol=1e-15;
[m,n]=size(r);
N=1:n;
F=find(x0>l+tol);
A=setdiff(N,F);
fx=zeors(m,1);
bx=zeors(m,1);
fx(F)=p(F);
bx(A)=p(A);
bx(bx>0)=0;
f1x=(x0-l)./a;
fxc=find(flx>fx);
f1x(fxc)=fx(fxc);
bxn=bx'*bx;
fxn=fx'*fx;
while bxn+fxn>delta
    if bxn<=L*f1x'*fx
        
    end
end