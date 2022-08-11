%Computes approximate solution of f(x)=0
%Input:function handle f; a,b such that f(a)*f(b)<0,
%        and tolerance tol
%Output: Approximate solution xc
function xc=bisect(f,a,b,tol,maxIt)
fa = f(a);
fb = f(b);
iter = 0;
if abs(fa) < tol
    xc = a;
    return;
end
if abs(fb)< tol
    xc = b;
    return;
end
if sign(fa)*sign(fb) >= 0 || a ==inf || b == inf
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
while (b-a)/2>tol && iter <maxIt
    iter = iter +1;
    c=(a+b)/2;
    fc=f(c);
    if fc == 0              %c is a solution, done
        break
    end
    if sign(fc)*sign(fa)<0  %a and c make the new interval
        b=c;
    else%c and b make the new interval
        a=c;
        fa=fc;
    end
end
xc=(a+b)/2;

%nthroot((0.5+sqrt(0.25+1/27)),3)+nthroot((0.5-sqrt(0.25+1/27)),3)