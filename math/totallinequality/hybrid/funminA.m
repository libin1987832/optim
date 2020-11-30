function f = funminA(x)
xn = x.x0 + x.alpha*x.u;
xn( xn < 0) = 0;
f = b - A * xn;
f(f<0) = 0;
f = sqrt(f'*f);