% computation A(rt>0,:)*p(n) A(rt>0,:)'*r(rt>0)(n)
function [steplength,retcode] = spiecewise(A,b,p,x,minrange,maxrange)
if p' * p < 1e-15
    steplength = 0;
    retcode = [2,0];
    return
end
retcode = [1,-1];
r = b - A * x;
Ap = A * p;
knot = r ./ Ap;
if min(knot) > maxrange || max(knot) < minrange
    knot = [minrange maxrange];
else
    knot = sort( knot( knot >= minrange & knot <= maxrange));
end
 for i = 1:length(knot)-1
     alpha = knot(i);
     beta = knot(i+1);
    % middle point give the active set
    rknot = r - 0.5 * (alpha+beta) * Ap;
    % Ap
    Ap = A( rknot > 0,:) * p;
    % the derive Ap'(r-alpha*Ap)=0
    Ar = A( rknot > 0 , : )' * r( rknot > 0);
    % alpha = (Ap*r)^T/(Ap'*Ap)
    steplength =( p' * Ar )/( Ap' * Ap );
    % the right point
    if steplength<=beta && steplength > alpha
        retcode = [1,0];
        break;
    else
        steplength = beta;
    end
end% if alph <-0.5
%     alph
% end
steplength = 0;
retcode = [3,0];