function alph=piecewise(A,b,p,x)
r=b-A*x;
ap=A*p;
ai=r./ap;
as=sort(ai(ai>0));
[am,an]=size(as);
for i =1:am
    t=as(i);
    if t<1
     rt=r-t*ap;
     rn=rt;
     rn(rt<0)=0;
     apn=ap;
     apn(rt<0)=0;
     alph=A'*rn/A'*apn;
        if alph<=t && alph > as(i-1)
            break;
        end
    else
        alph=1;
        break;
    end
end
