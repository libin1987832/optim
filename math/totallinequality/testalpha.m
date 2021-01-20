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
pxy(1).Y = knoty;烦烦烦
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

function xk=MPRGP(A,b,x0,yk,L,a,delta,maxIter)
tol=1e-15;
Ftol = 1e-10;
[m,n]=size(A);
g = A' * yk;
F = x0 > Ftol;
fx = zeros(n , 1);
bx = zeros(n , 1);
fx = g(F);
bx = min(g(~F),0);
iter = 0;
bk = yk + A*x0;
r = yk;
p = fx;
while bx' * bx + fx' * fx > delta && iter < maxIter
    iter = iter + 1;
    f1x = min( x0(F) / a , fx(F));
    if bxn <= L * f1x' * fx
        Ap = A*p; 
        acg = (r' * Ap) / (Ap' * Ap); y = x0 - acg * p;
        I = p > Ftol;
        af = min(x0(I) ./ p(I));
        if acg < af
            xk = y; 
            r = r - acg * Ap;
            g = A' * r;
            F = xk > Ftol;
            fx = 0;
            fx = g(F);
            Afx = A' * fx;
            grma=(Afx' * Ap)/(Ap' * Ap);
            p = fx - grma * p;
            disp(['conject gradient:',num2str(r'*r)]);
        else
            afp = af * p;
            xk2 = x0 - afp;
            r = r - A * afp;
            F = xk2 > Ftol;
            g = A' * r;
            fx = 0;
            fx = g(F);
            xk=xk2-a*fx;
            F = xk > Ftol;
            xk( ~F )=0;
            r = A * xk - bk;
            g = A' * r;
            fx = 0;
            fx = g(F);
            p = fx;
            disp(['out range in the subspace:',num2str(r'*r)]);
        end
    else
        d=bx;
        acg=(r'*d)/(d'*A*d);
        xk=x0-acg*d;r=r-acg*A*d;[p,bx]=fbxf(xk,m,l,tol,r);
        fx=p;
        disp(['go to implement space:',num2str(r'*r)]);
    end
    k=k+1;
    bxn=bx'*bx;
    x0=xk;
end
end

