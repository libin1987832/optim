% computation A(rt>0,:)*p(n) A(rt>0,:)'*r(rt>0)(n)
function alph=spiecewise(A,b,p,x,range)
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
    % if t is greater than one, exit
    if t<range
    % middle point give the active set    
     rt=r-0.5*(t+tas(i))*ap;
     % Ap
     Ad=A(rt>0,:)*p;
     % the derive Ap'(r-alpha*Ap)=0
     Ar=A(rt>0,:)'*r(rt>0);
     % alpha = (Ap*r)^T/(Ap'*Ap)
     alph=(p'*Ar)/(Ad'*Ad);
     % the left point
     if i>2
        last=as(i-1);
     else
       last=0;
     end
     % the right point
        if alph<=t && alph > last
            break;
        end
    else
        alph=range;
        break;
    end    
end
% if alph <-0.5
%     alph
% end
