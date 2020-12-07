% computation A(rt>0,:)*p(n) A(rt>0,:)'*r(rt>0)(n)
function [steplength, retcode] = spiecewise(A,b,p,x,minrange,maxrange)
display = true;
if p' * p < 1e-15
    steplength = 0;
    retcode = [2,0];
    return
end
r = b - A * x;
Ap = A * p;
knot = r ./ Ap;
knot = sort( knot( knot >= minrange & knot <= maxrange));
knot = unique([minrange;knot;maxrange]);
for i = 1:length(knot)-1
    alpha = knot(i);
    beta = knot(i+1);
    % middle point give the active set
    rknot = r - 0.5 * (alpha + beta) * Ap;
    % Ap
    Apr = Ap(rknot > 0);
    % the derive Ap'(r-alpha*Ap)=0
    Ar = Apr' * r(rknot > 0);
    % alpha = (Ap*r)^T/(Ap'*Ap)
    steplength = Ar / ( Apr' * Apr );
    % the right point
    if steplength >= beta
        steplength = beta;
        retcode = [2,1];
    elseif steplength <= alpha
        steplength = alpha;
        retcode = [1,1];
        break;
    else
        retcode = [1,0];
        break;
    end
end
% if ~exist('steplength','var')
%     t=1;
% end