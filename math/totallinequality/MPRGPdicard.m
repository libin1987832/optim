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

