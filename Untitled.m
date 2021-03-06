function test()
m=10;n=10;
A = rand(m,n);
b = rand(m,1);
p = rand(n,1);
x0 = rand(n,1);
xa = [-10:0.01:10];
ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
 pxy={};
% pxy(1).X = knot;
% pxy(1).Y = knoty;
pxy(1).X = xa;
pxy(1).Y = ya;
figure
hold on
p1 = arrayfun(@(a) plot(a.X,a.Y),pxy,'UniformOutput',false);


hold off
 save('nonconvexalpha.mat','A','b','p','x0','xa','ya')
end
function f = funmin(A,b,x0,u,alpha)
xn = x0 + alpha*u;
xn( xn < 0) = 0;
f = b - A * xn;
f(f<0) = 0;
f = 0.5*(f'*f);
end