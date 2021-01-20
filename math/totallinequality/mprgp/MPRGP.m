function xk=MPRGP(A, b, x0, L, a, delta, Ftol, maxIter)
[~,n]=size(A);
Ax0 = A*x0;
yk = b - Ax0;
yk(yk<0) = 0;
g = A' * yk;
F = x0 > Ftol;
fx = zeros(n , 1);
bx = zeros(n , 1);
f1x = zeros(n , 1);
fx( F ) = g( F );
bx( ~F ) = min(g( ~F ),0);
iter = 0;
bk = yk + Ax0;
r = yk;
p = fx;
while bx' * bx + fx' * fx > delta && iter < maxIter
    iter = iter + 1;
    f1x( : ) = 0;
    f1x( F ) = min( x0( F ) / a , fx( F ));
    if bx' * bx <= L * f1x' * fx
        Ap = A * p;
        acg = (r' * Ap) / (Ap' * Ap); y = x0 - acg * p;
        I = p > Ftol;
        af = min(x0(I) ./ p(I));
        if acg < af
            xk = y;
            r = r - acg * Ap;
            g = A' * r;
            F = xk > Ftol;
            fx( : ) = 0;
            fx( F ) = g(F);
            Afx = A' * fx;
            grma=(Afx' * Ap)/(Ap' * Ap);
            p = fx - grma * p;
            disp(['conject gradient:',num2str(r'*r)]);
        else
            afp = af * p;
            xk2 = x0 - afp;
            r = A * xk2 - bk;
            F = xk2 > Ftol;
            g = A' * r;
            fx( : ) = 0;
            fx( F ) = g(F);
            xk = xk2 - a * fx;
            F = xk > Ftol;
            xk( ~F )=0;
            r = A * xk - bk;
            g = A' * r;
            fx( : ) = 0;
            fx( F ) = g( F );
            p = fx;
            disp(['out range in the subspace:',num2str(r'*r)]);
        end
    else
        Ad = A * bx;
        acg = (Ad' * r) / (Ad' * Ad);
        xk = x0 - acg * d; r = A * xk - bk;
        g = A' * r;
        F = xk > Ftol;
        p = g(F);
        fx = p;
        disp(['go to implement space:',num2str(r'*r)]);
    end
    bx( : ) = 0;
    bx( ~F ) = min( g( ~F ) , 0);
    x0 = xk;
end
end

