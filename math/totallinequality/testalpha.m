addpath('./hybrid');
A = [1 -2];
u = [-1;-2];
x0 = [1;1];
b = [10];
[alpha, minf, knot] = arraySpiecewise(A,b,x0,u);
[alpha1, minf, knot1] = arraySpiece(A,b,x0,u);
%assert(~sum(knot - [0,1/14,3/16,5/18,1/3,1/3+3/27,1/3+7/24,1/3+11/21,1]));
knoty = arrayfun(@(alpha) funmin(A,b,x0,u,alpha), knot);
xa = [0:0.01:1];
ya = arrayfun(@(alpha) funmin(A,b,x0,u,alpha), xa);
pxy={};
pxy(1).X = knot;
pxy(1).Y = knoty;
pxy(2).X = xa;
pxy(2).Y = ya;
figure
hold on
p = arrayfun(@(a) plot(a.X,a.Y),pxy);
p(1).Marker = 'o';
p(2).Marker = '+';
hold off
alpha, minf,
function f = funmin(A,b,x0,u,alpha)
xn = x0 + alpha*u;
xn( xn < 0) = 0;
f = b - A * xn;
f(f<0) = 0;
f = 0.5*(f'*f);
end