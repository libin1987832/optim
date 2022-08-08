%Computes approximate solution of f(x)=0
%Input:function handle f; a,b such that f(a)*f(b)<0,
%        and tolerance tol
%Output: Approximate solution xc
function xc=bisect(f,a,b,epsbisect,fa,fb)
if abs(fa) < epsbisect
    xc = a;
    return;
end
if abs(fb)< epsbisect
    xc = b;
    return;
end
if sign(fa)*sign(fb) >= 0
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
while (b-a)/2>epsbisect
    xc=(a+b)/2;
    fc=f(xc);
    if fc == 0              %c is a solution, done
        break
    end
    if sign(fc)*sign(fa)<0  %a and c make the new interval
        b=xc;
    else                  %c and b make the new interval
        a=xc;
        fa=fc;
    end
end
end

