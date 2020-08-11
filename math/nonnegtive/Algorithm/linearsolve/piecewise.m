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
     Ad=A(rt>0,:)*p;
     Ar=A(rt>0,:)'*r(rt>0);
     alph=(p'*Ar)/(Ad'*Ad);
     if i>2
        last=as(i-1);
     else
       last=0;
     end
        if alph<=t && alph > last
            break;
        end
    else
        alph=1;
        break;
    end
end
