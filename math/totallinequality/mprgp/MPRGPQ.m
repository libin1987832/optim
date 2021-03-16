% x^TAx+bx
%x(k+1) = x(k) - ad(gradient g); ||Ax - bk ||
function [x0,rpk] = MPRGPQ(A, bk, x0, L_cons, a, delta, Ftol, maxIter,bnorm)
[~,n]=size(A);
debug = 2;%1 precondition %2 debug
r = A*x0-bk;
F = x0 > Ftol;
fx = zeros(n , 1);
bx = zeros(n , 1);
f1x = zeros(n , 1);
bx( ~F ) = min(r( ~F ),0);
iter = 0;
fx(F)=r(F);
p = fx;
update = true;
% if debug
%
%     bnorm=bk'*bk;
% end
while bx' * bx + fx' * fx > delta && iter < maxIter
    iter = iter + 1;
    f1x( : ) = 0;
    f1x( F ) = min( x0( F ) / a , fx( F ));
    if bx' * bx <= L_cons * f1x' * fx
        if debug == 2 && update
            AF=A(F,F);
            nf=size(F,1);
            L=ichol(AF,struct('michol','on'));
            y=L'\(L\r(F));
            p(F)=y;
        end
        Ap = A * p;
        rp=r'*p;
        acg = rp / (p' * Ap); y = x0 - acg * p;
        I = p > Ftol;
        if any(I)
            af = min(x0(I) ./ p(I));
        else
            af = +inf;
        end
        if acg < af
            %iter = iter - 1;
            xk = y;
            r = r - acg * Ap;
            fx( F ) = r( F );
            if debug == 2
                AF=A(F,F);
                L=ichol(AF,struct('michol','on'));
                y=L'\(L\r(F));
                grma=(r(F)'*y)/rp;
            else
                grma=(fx' * Ap) / (p' * Ap);
            end
            %             F = xk > Ftol;
            %             fx( : ) = 0;
            %             fx( F ) = r(F);
            
            if debug
                r_d0 = x0'*A*x0-2*x0'*bk+bnorm;
                r_0 = xk'*A*xk-2*xk'*bk+bnorm;
                disp(['old_r ',num2str(r_d0),' conject gradient:',num2str(r_0),' af: ', num2str(af), ' acg:',num2str(acg)]);
            end
            p(F) = y - grma * p(F);
            update = false;
        else
            p = fx;
            Ap=A*p;
            afp = af * p;
            xk2 = x0 - afp;
            rold = r;
            r = r - af*Ap;
            F = xk2 > Ftol;
            fx( : ) = 0;
            fx( F ) = r(F);
            xk = xk2 - a * fx;
            %             xk = x0 - a * p;
            F = xk > Ftol;
            xk( ~F )=0;
            r = A * xk - bk;
            fx( : ) = 0;
            fx( F ) = r( F );
            if debug
                r_d0 = x0'*A*x0-2*x0'*bk+bnorm;
                r_0 = xk'*A*xk-2*xk'*bk+bnorm;
                disp(['old_r ',num2str(r_d0),' projected:',num2str(r_0),' af: ', num2str(af), ' acg:',num2str(acg),' a:',num2str(a),' p*f>0:',num2str(rold'*p)]);
            end
            p = fx;
            update = true;

            %              break;
        end
    else
        d=bx;
        Ad = A * bx;
        acg = (d' * r) / (d' * Ad);
        xk = x0 - acg * bx; r = r-acg * Ad;
        F = xk > Ftol;
        fx( : ) = 0;
        fx( F ) = r( F );
        p = fx;
        if debug
            r_d0 = x0'*A*x0-2*x0'*bk+bnorm;
            r_0 = xk'*A*xk-2*xk'*bk+bnorm;
            disp(['old_r ',num2str(r_d0),' implement:',num2str(r_0), ' acg:',num2str(acg)]);
        end
        update = true;
        %          break;
    end
    bx( : ) = 0;
    bx( ~F ) = min( r( ~F ) , 0);
    x0 = xk;
end
rpk = bk-A*x0;
end
