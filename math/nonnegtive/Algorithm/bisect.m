function c=bisect(r,p)
TOL = 1e-10;
a=0;
b=1;
while f(r,p,b)<0
    b=b/2;
end
while (b-a)/2>TOL
    c=(a+b)/2;
    if f(r,p,c) == 0
        break;
    end
    if f(r,p,a)*f(r,p,c)<0
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
function [am,bm]=sec(r,p,a,b)
c=(a+b)/2;
if f(r,p,c)>0
    am=a;
    bm=c;
    return;
else
    [am,bm]=sec(a,c);
    [am,bm]=sec(c,b);
    return;
end
end