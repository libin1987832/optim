function test()
addpath('./hybrid')
p=0;
load('test.mat','A','b','x0');
while  norm(p) < 1e-10
m=20;n=10;
% A = rand(m,n)-0.5;
% b = rand(m,1);
% x = rand(n,1);
% x0 = rand(n,1);
% A = A;
% b = b.*10;
% x = -x.*10;
% b = A*x-0.001;
r = (b-A*x0);
r(r<0) = 0;
p = (A'*A)\(A'*r);
%p = A'*r;
end
[x0,p,x0+0.0530*p,x0+0.0540*p -x0./p]
x1 = x0+0.0530*p;
x1(x1<0) = 0;
f1 = A*(x1-x0);
f1 = norm(f1);
x2 = x0+0.0540*p;
x2(x2<0) = 0;
f2 = A*(x2-x0);
f2 = norm(f2);

[alpha, minf, knot,retcode] = arraySpiece(A,b,x0,p,1e-10,300);
[f1 f2 f1/0.0530 f2/0.0540 alpha funmin(A,b,x0,p,alpha-0.001),funmin(A,b,x0,p,alpha),funmin(A,b,x0,p,alpha+0.001),]
xf1 = (alpha-0.001)*p+x0;
If1 = xf1<0; 
xf1(If1) = 0;
xf2 = (alpha+0.001)*p+x0;
If2 = xf2<0; 
xf2(If2) = 0;
ff1 = b-A*xf1;
ff1(ff1<0)=0;
Aff1 = A'*ff1;
ff2 = b-A*xf2;
ff2(ff2<0)=0;
Aff2 = A'*ff2;
p1 = p;
p1(If1) = 0;
p2 = p;
p2(If2) = 0;

[-Aff1'*p1 -Aff2'*p2]
%  save('test.mat','A','b','x0');
f0 = funmin(A,b,x0,p,0);
xa = [0.3:0.001:0.35];
ya = arrayfun(@(alpha)  funmin(A,b,x0,p,alpha), xa);

ya2 = arrayfun(@(alpha) funnorm(A,x0,p,alpha), xa);
dn =60;
% [ya2(1:dn);
% ya(1:dn);
% xa(1:dn)]

mya = mean(ya);
ya(abs(ya-mya)<1e-10) = mya;
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
f =  A*(xn-x0);
f = sqrt((f'*f))/alpha;
end

