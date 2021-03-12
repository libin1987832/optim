% x^TAx+bx
%x(k+1) = x(k) - ad(gradient g); ||Ax - bk ||
function [x0,rpk] = MPRGPQ(A, bk, x0, L, a, delta, Ftol, maxIter)
[~,n]=size(A);

r = A*x0-bk;
F = x0 > Ftol;
fx = zeros(n , 1);
bx = zeros(n , 1);
f1x = zeros(n , 1);
fx( F ) = r( F );
bx( ~F ) = min(r( ~F ),0);
iter = 0;
p = fx;
debug = 0;
while bx' * bx + fx' * fx > delta && iter < maxIter
    iter = iter + 1;
    f1x( : ) = 0;
    f1x( F ) = min( x0( F ) / a , fx( F ));
    if bx' * bx <= L * f1x' * fx
        Ap=A*p;
        acg = (r' * p) / (p' * Ap); y = x0 - acg * p;
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
            F = xk > Ftol;
            fx( : ) = 0;
            fx( F ) = r(F);
            grma=(fx' * Ap) / (p' * Ap);
            if debug
                p1 = fx - grma * p;
                r_d0 = A*x0-bk;
            disp(['old_r ',num2str(0.5*r_d0'*r_d0),' conject gradient:',num2str(0.5*r'*r),' diff: ', num2str(0.5*r_d0'*r_d0-0.5*r'*r),' p1*A*Ap:',...
                num2str(p1'*A'*A*p),' cos(p1,p)',...
                num2str(p1'*p/(norm(p1)*norm(p))),' af: ', num2str(af),...
                ' acg:',num2str(acg)]);
            end
            p = fx - grma * p;
        else
            Ap=A*p;
            afp = af * p;
            xk2 = x0 - afp;
            r = r - af*Ap;
            F = xk2 > Ftol;
            fx( : ) = 0;
            fx( F ) = r(F);
            xk = xk2 - a * fx;
            if debug
                r_d0 = A*x0-bk;
                g_d = A'*(A*x0-bk);
                disp(['old: ',num2str(0.5*r_d0'*r_d0),' new: ',num2str(0.5*r'*r),'p*g:',num2str(p'*g_d)]);            
            end
%             xk = x0 - a * p;
            F = xk > Ftol;
            xk( ~F )=0;
            r = A * xk - bk;
            fx( : ) = 0;
            fx( F ) = r( F );
            p = fx;
            if debug
  %              x2 = pcg(A,b,[],100,[],[],x0);
%               xk_d = x0 - 0.000001 * p;r_d2 = A*xk_d-bk;
%               disp(['r_d:',num2str(r_d2'*r_d2)]);
                g_d = A'*(A*x0-bk);
              disp([' af: ', num2str(af),...
                ' acg:',num2str(acg),'old: ', num2str(0.5*r_d0'*r_d0), ...
                ' out range in the subspace:',num2str(0.5*r'*r ),'diff: ',num2str(0.5*r_d0'*r_d0-0.5*r'*r),' xa-x,fx', num2str(-(xk-x0)'*g_d)]);
            end
%              break;
        end
    else
        Ad = A * bx;
        acg = (d' * r) / (d' * Ad);
        xk = x0 - acg * bx; r = r-acg * Ad;
        F = xk > Ftol;
        fx( : ) = 0;
        fx( F ) = r( F );
        p = fx;
        if debug
          r_d0 = A * x0 - bk;
          disp(['old ',num2str(0.5*r_d0'*r_d0),' go to implement space:',num2str(0.5*r'*r)]);
        end
%          break;
    end
    bx( : ) = 0;
    bx( ~F ) = min( r( ~F ) , 0);
    x0 = xk;
end
rpk = bk-A*x0;
end
