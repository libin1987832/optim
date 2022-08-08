function [d, lambda, dAx, Md] = searchdir(M, Ax, xB, tao, x, method, tol)
    xBx = xB * x;
    if method == 1
       options = optimoptions('quadprog','Display','off');%为了不输出 quadprog中间的过程
       [d,~,~,~,~]=quadprog(M, Ax, [], [], xB, tao, -x, [],[],options);
       Md = M * d;
       xMd = x' * Md;
       dMd = d' * Md;
       dAx = d' * Ax;
       dBx = xB * d;
       xAx = x' * Ax;
       lambda = ( xMd + xAx + dMd + dAx )/( xBx + dBx);
    elseif method == 2
         Bx = xB';
         Mx = M .* x;
         invM = 1 ./ M;
         Bx2 = Bx .^ 2;
         
         
         
         grt0Bx = Bx > 0;
         Lam0=( Ax( grt0Bx ) - Mx( grt0Bx ) ) ./ Bx( grt0Bx );
         lambdaLower = min( Lam0 );
         
         Lam1 =  ( ( invM .* Bx )' * Ax +tao ) / ( invM' * Bx2 );
         funDir = @(lam) Bx' * directionByD(lam, x, Bx, Ax, Mx, invM) -tao; 
         
         if funDir(Lam1)* funDir( lambdaLower)<1e-12
             lambdaUpper = Lam1; 
         else
         Lam2 = max( Lam0);
         grt1Bx = abs(Bx)> 0;
         Lam3 =( invM' * ( abs( Bx ) .* abs( Ax ) ) +tao) /min( invM(grt1Bx) .* Bx2(grt1Bx) );
         lambdaUpper = max([Lam1, Lam2, Lam3]);
         end
         
        lambda=bisect(funDir, lambdaLower, lambdaUpper, tol);
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