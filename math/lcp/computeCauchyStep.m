function [cauthy,xc]=computeCauchyStep(c,G,x0,d)
grad=G*x0+c;
dgradd=d'*G*d;
aglobal=-(d'*grad)/(dgradd);
assert(aglobal>0, ['aglobal = ' num2str(aglobal) '<0 is impossible!'])
grt0=(d<0);
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
%  test_cauchy(cauthy,d,G,c,x0)