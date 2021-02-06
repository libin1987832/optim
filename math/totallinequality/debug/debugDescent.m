function debugDescent(A,b,x0,p)
xa = [0:1e-2:5];
ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
pxy={};
pxy(1).X = xa;
pxy(1).Y = ya;
figure
hold on
p1 = arrayfun(@(a) plot(a.X,a.Y),pxy,'UniformOutput',false);
end
function f = funmin(A,b,x0,u,alpha)
xn = x0 + alpha*u;
xn( xn < 0) = 0;
f = b - A * xn;
f(f<0) = 0;
f = 0.5*(f'*f);
end