function [d, lambda] = searchdir(M, Ax, xB, tao, x, method, tol)
    xBx = xB' * x;
    if method == 1
       [d,fval,exitflag,output,lambda]=quadprog(M, Ax, [], [], xB, tao, -x, []);
       Md = M * d;
       xMd = x' * Md;
       dMd = d * Md;
       dAx = d' * Ax;
       dBx = d' * xB';
       lambda = ( xMd + xAx + dMd + dAx )/( xBx + dBx);
    elseif method == 2
         D = M;
         Bx = xB';
         Dx = D .* x;
         invD = 1 ./ D;
         c = 0.5 * xBx -1;
         Bx2 = Bx .^ 2;
         Lam1 =  ( invD .* Bx * Ax - c ) / ( invD' * Bx2 );  
         Lam2 = max((Ax - Dx) ./ Bx);
         Lam3 = invD' * ( abs( Bx ) .* abs( Ax ) ) - c /min( invD .* Bx2 ); 
         lambdaUpper = max(Lam1, Lam2, Lam3);
         grt0Bx = Bx > 0;
         lambdaLower = min( Ax( grt0Bx ) - Dx( grt0Bx ) ./ Bx( grt0Bx ));
         funDir = @(lam) Bx' * directionByD(lam, Bx, Ax, Dx) + c;
         lambda=bisect(funDir, lambdaUpper, lambdaLower, tol);
         d = directionByD(lambda, Bx, Ax, Dx);
    end
end
function dir = directionByD(lambda, Bx, Ax, Dx)
         lBx = lambda * Bx;
         lBxAx = lBx - Ax;
         grt0lBAx = lBxAx > Dx; 
         dir = -x;
         dir( grt0lBAx ) = lBxAx( grt0lBAx );
end