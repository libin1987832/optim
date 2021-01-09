function [A,b,x0]=test()
p=0;
while  norm(p) < 1e-10
m=2;n=2;
A = rand(m,n);
b = rand(m,1);
x0 = rand(n,1);
A = A;
b = b.*10;
x0 = x0.*10;
r = (b-A*x0);
r(r<0) = 0;
p = (A'*A)\(A'*r);
p = A'*r;
end
f0 = funmin(A,b,x0,p,0);
xa = [0.0001:0.1:2];
ya2 = arrayfun(@(alpha)  funmin(A,b,x0,p,alpha), xa);
ya2(1:5)
ya = arrayfun(@(alpha) funnorm(A,x0,p,alpha), xa);
ya(10:15)
 pxy={};
% pxy(1).X = knot;
% pxy(1).Y = knoty;
pxy(1).X = xa;
pxy(1).Y = ya;
% pxy(2).X = xa;
% pxy(2).Y = ya2;
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
function f = funnorm(A,x0,u,alpha)
xn = x0 + alpha*u;
xn( xn < 0) = 0;
f =  (xn-x0);
f = sqrt((f'*f))/alpha;
end
save('test.mat','A','b','x0');