function [d, lambda, dAx, Md] = searchdir(M, Ax, xB, tao, x, method, tol)
    xBx = xB * x;
    if method == 1
        opts = optimoptions('quadprog','Display','off');
       d=quadprog(M, Ax, [], [], xB, tao, -x, [], [], opts);
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
         c = 0.5 * xBx -1;
         Bx2 = Bx .^ 2;
         Lam1 =  ( ( invM .* Bx )' * Ax - c ) / ( invM' * Bx2 );  
%          Lam2 = max((Ax - Mx) ./ Bx);
%          Lam3 = (invM' * ( abs( Bx ) .* abs( Ax ) ) - c )/min( invM .* Bx2 ); 
%          lambdaUpper = max([Lam1, Lam2, Lam3]);
        lambdaUpper = Lam1;
         grt0Bx = Bx > 0;
         lambdaLower = min( ( Ax( grt0Bx ) - Mx( grt0Bx ) ) ./ Bx( grt0Bx ));
         funDir = @(lam) Bx' * directionByD(lam, x, Bx, Ax, Mx, invM) + c;
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