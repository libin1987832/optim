function test()
addpath('./hybrid')
p=0;
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
load('test2.mat','A','b','x0');
r = (b-A*x0);
r(r<0) = 0;
p = (A'*A)\(A'*r);
%p = A'*r;
end
% save('test2.mat','A','b','x0');
% [x0,p,x0+0.0530*p,x0+0.0540*p -x0./p]
aa = -x0./p;
aa = sort(aa>=0)
% test some point ||A(x(a)-x)||
x1 = x0+0.0530*p;
x1(x1<0) = 0;
f1 = A*(x1-x0);
f1 = norm(f1);
x2 = x0+0.0540*p;
x2(x2<0) = 0;
f2 = A*(x2-x0);
f2 = norm(f2);

[alpha, minf, knot,retcode] = arraySpiece(A,b,x0,p,1e-10,300);
alpha

% compute the gradient of function || (b-Ax)_+ ||
xf0 = alpha*p+x0;
If0 = xf0 < 0;
xf0(If0) = 0;
ff0 = b-A*xf0;
ff0(ff0<0)=0;
Aff0 = A'*ff0;
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
p0 = p;
p0(If0) = 0;
[-Aff1'*p1 -Aff2'*p2  -Aff0'*p0]

f0 = funmin(A,b,x0,p,0);
xa = [0.001:0.001:1.5];
ya = arrayfun(@(alpha)  funmin(A,b,x0,p,alpha), xa);

ya2 = arrayfun(@(alpha) f0 + 0.5*funnorm(A,x0,p,alpha), xa);
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
pxy(2).X = xa;
pxy(2).Y = ya2;
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

