addpath('./hybrid');
A = [1 4 7;2 5 8;3 6 9];
u = [-1;2;-3];
x0 = [1;1;1];
b = [11;12;13];
[alpha, minf, knot] = arraySpiecewise(A,b,x0,u);
assert(~sum(knot - [0,1/14,3/16,5/18,1/3,1/3+3/27,1/3+7/24,1/3+11/21,1]));
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