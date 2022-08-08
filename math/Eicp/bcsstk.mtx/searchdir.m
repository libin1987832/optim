function [d, lambda, dAx, Md] = searchdir(M, Ax, xB, tao, x, method, epsbisect)
   
    if method == 1
       opts = optimoptions('quadprog','Display','off');
       d=quadprog(M, Ax, [], [], xB, tao, -x, [], [], opts);
       Md = M * d;
       xMd = x' * Md;
       dMd = d' * Md;
       dAx = d' * Ax;
       dBx = xB * d;
       xAx = x' * Ax;
       lambda = ( xMd + xAx + dMd + dAx )/( xB * x + dBx);
    elseif method == 2
         Bx = xB';
         Mx = M .* x;
         invM = 1 ./ M;
         Bx2 = Bx .^ 2;
         funDir = @(lam) Bx' * directionByD(lam, x, Bx, Ax, Mx, invM) -tao; 
         
         grt0Bx = Bx > 0;
         Lam0=( Ax( grt0Bx ) - Mx( grt0Bx ) ) ./ Bx( grt0Bx );
         lambdaLower = min( Lam0 );
         fa=funDir( lambdaLower);
         
         Lam1 =  ( ( invM .* Bx )' * Ax +tao ) / ( invM' * Bx2 );
         fb=funDir(Lam1);
         
         if fa* fb<1e-10
             lambdaUpper = Lam1; 
         else
         Lam2 = max( Lam0);
         grt1Bx = abs(Bx)> 0;
         Lam3 =( invM' * ( abs( Bx ) .* abs( Ax ) ) +tao) /min( invM(grt1Bx) .* Bx2(grt1Bx) );
         lambdaUpper = max([Lam1, Lam2, Lam3]);
         fb=funDir( lambdaUpper);
         end
         
        lambda=bisect(funDir, lambdaLower, lambdaUpper, epsbisect,fa,fb);
        d = directionByD(lambda, x, Bx, Ax, Mx, invM);


         
         dAx = d' * Ax;
         Md = M .* d;
    end
end
function dir = directionByD(lambda, x, Bx, Ax, Dx, invM)
         lBx = lambda * Bx;
         lBxAx = lBx - Ax;
         grt0lBAx = lBxAx > -Dx; 
         dir = -x;
         dir( grt0lBAx ) = invM( grt0lBAx ) .* lBxAx( grt0lBAx );
end