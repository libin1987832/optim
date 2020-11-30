% computation A(rt>0,:)*p(n) A(rt>0,:)'*r(rt>0)(n)
function steplength = spiecewise(A,b,p,x,alpha,beta)
if p' * p < 1e-15
    steplength = 0;
    return
end
r = b - A * x;
Ap = A * p;
knot = r ./ Ap;
knot = sort( knot( knot > 0 ));
% middle point give the active set
rknot = r - 0.5 * (alpha+beta) * Ap;
% Ap
Ap = A( rknot > 0,:) * p;
% the derive Ap'(r-alpha*Ap)=0
Ar = A( rknot > 0 , : )' * r( rknot > 0);
% alpha = (Ap*r)^T/(Ap'*Ap)
alph=( Ap' * Ar )/( Ap' * Ap );
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
end% if alph <-0.5
%     alph
% end
