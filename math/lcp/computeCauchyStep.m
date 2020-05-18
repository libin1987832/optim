function [cauthy,xc,res]=computeCauchyStep(c,G,x0,d)
tol=1e-10;
if norm(d)<tol
    res=0;
    xc=x0;
    cauthy=0;
    return;
end
grad=G*x0+c;
dgradd=d'*G*d;
aglobal=-(d'*grad)/(dgradd);
assert(aglobal>0, ['aglobal = ' num2str(aglobal) '<0 is impossible!'])
% mark less 0 index
grt0=(d<0);
% if all indredient great 0 then natural sartisfied
if sum(grt0)<1
    amax=aglobal
else
    amax=min(x0(grt0)./-d(grt0));
end
if aglobal<amax 
    cauthy=aglobal;
else
    cauthy=amax;
end
 xc=x0+cauthy*d;
res=1;
%  test_cauchy(cauthy,d,G,c,x0)