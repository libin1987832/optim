% || Ax - yk - Axk || || aAd - (-yk) || ||aAd - rk ||
%x(k+1) = x(k) - ad(gradient g); ||Ax - bk ||
function [x0,rpk] = MPRGPlsqr(A, b, x0,rpk, L, a, delta, Ftol, maxIter)
[~,n]=size(A);
% Ax0 = A*x0;
% yk = b - Ax0;
yk = rpk;
zk = -rpk;
yk(yk<0) = 0;
zk(zk<0) = 0;
% bk = yk + Ax0;
% r = Ax0 - bk;
bk = zk + b;
r = -yk;
g = A' * r;
F = x0 > Ftol;
fx = zeros(n , 1);
bx = zeros(n , 1);
f1x = zeros(n , 1);
fx( F ) = g( F );
bx( ~F ) = min(g( ~F ),0);
iter = 0;
p = fx;
debug = 1;
while bx' * bx + fx' * fx > delta && iter < maxIter
    iter = iter + 1;
    f1x( : ) = 0;
    f1x( F ) = min( x0( F ) / a , fx( F ));
    if bx' * bx <= L * f1x' * fx
        Ap = A * p;
        acg = (r' * Ap) / (Ap' * Ap); y = x0 - acg * p;
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
            g = A' * r;
            F = xk > Ftol;
            fx( : ) = 0;
            fx( F ) = g(F);
            Afx = A * fx;
            grma=(Afx' * Ap) / (Ap' * Ap);
            if debug
                 AF = A(:,F);
                 xkkr = xk;
%                 L = ichol(AF'*AF,struct('michol','on'));
%                 [x2,fl2,rr2,it2,rv2] = pcg(AF'*AF,AF'*bk,1e-8,3,L,L');
                bkkr = yk + A*x0; 
                xkkr(F) = krylovk(AF,bkkr,3);
                [rpkkr, normrkr, xminkr, Ar, KKTkr, face1, face2] = kktResidual(A, b, xkkr , [],1);
                p1 = fx - grma * p;
                r_d0 = A*x0-bk;
            disp(['old_r ',num2str(0.5*r_d0'*r_d0),' conject gradient:',num2str(0.5*r'*r),' diff: ', num2str(0.5*r_d0'*r_d0-0.5*r'*r),' p1*A*Ap:',...
                num2str(p1'*A'*A*p),' cos(p1,p)',...
                num2str(p1'*p/(norm(p1)*norm(p))),' af: ', num2str(af),...
                ' acg:',num2str(acg),'KKTkr: ',num2str(KKTkr)]);
            end
            p = fx - grma * p;
        else
            afp = af * p;
            xk2 = x0 - afp;
            r = A * xk2 - bk;
            F = xk2 > Ftol;
            g = A' * r;
            fx( : ) = 0;
            fx( F ) = g(F);
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
            g = A' * r;
            fx( : ) = 0;
            fx( F ) = g( F );
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
        acg = (Ad' * r) / (Ad' * Ad);
        xk = x0 - acg * bx; r = A * xk - bk;
        g = A' * r;
        F = xk > Ftol;
        fx( : ) = 0;
        fx( F ) = g( F );
        p = fx;
        if debug
          r_d0 = A * x0 - bk;
          disp(['old ',num2str(0.5*r_d0'*r_d0),' go to implement space:',num2str(0.5*r'*r)]);
        end
%          break;
    end
    bx( : ) = 0;
    bx( ~F ) = min( g( ~F ) , 0);
    x0 = xk;
end
rpk = -r - zk;


end
