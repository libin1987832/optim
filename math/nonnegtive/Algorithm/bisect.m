function c=bisect(r,p)
TOL = 1e-5;
a=0;
b=1;
r(r<0)=0;
f0=p'*r;
r1=r-p;
r1(r1<0)=0;
f1=p'*r1;
if abs(f1)<1e-10 || f1>0
    c=1;
    return;
end
% [b,fs,s]=sec(r,p,0,1);
% if abs(fs)<1e-5
%     c=b;
%     return;
% end
while (b-a)/2>TOL
    c=(a+b)/2;
    v=f(r,p,c);
    if  abs(v)<1e-5
        break;
    end
    va=f(r,p,a);
    vb=f(r,p,c);
    if vb*va<0
     b=c;
    else
        a=c;
    end
end
end
function v=f(r,p,a)
    rp=r-a*p;
    rp(rp<0)=0;
    v=p'*rp;
end
function [bm,fk,s]=sec(r,p,a,b)
c=(a+b)/2;
fk=f(r,p,c);
if fk>0
    bm=c;
    s=0;
elseif  abs(a-b)<1e-5
    bm=1;
    s=1;
else
    [bm,fk,s]=sec(r,p,c,b);
    if s==1
        [bm,fk,s]=sec(r,p,a,c);
    end
end
end