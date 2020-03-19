function c=bisect(A,B,p,x)
TOL = 1e-10;
a=0;
b=1;
while f(A,B,p,x,b)>0
    b=b/2;
end
while (b-a)/2>TOL
    c=(a+b)/2;
    if f(A,B,p,x,c) == 0
        break;
    end
    if f(A,B,p,x,a)*f(A,B,p,x,c)<0
     b=c;
    else
        a=c;
    end
end
end
function v=f(A,B,p,x,a)
    rp=B-A*x-a*A*p;
    rp(rp<0)=0;
    v=2*p'*A'*rp;
end