% computation A(rt>0,:)*p A(rt>0,:)'*r(rt>0)
function alph=spiecewise(A,b,p,x)
r=b-A*x;
if p'*p<1e-15
    alph = 0;
    return
end
ap=A*p;
ai=r./ap;
as=sort(ai(ai>0));
tas=[0;as];
[am,an]=size(as);
%alph=-1;
for i =1:am
    t=as(i);
    if t<1
     rt=r-0.5*(t+tas(i))*ap;
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
% if alph <-0.5
%     alph
% end
