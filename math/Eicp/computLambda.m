function lambda = computLambda(M, Ax, xAx, xB, xBx, d, x, method)
  if method == 1
       Md = M * d;
       xMd = x' * Md;
       dMd = d * Md;
       dAx = d' * Ax;
       dBx = d' * xB';
       lambda = ( xMd + xAx + dMd + dAx )/( xBx + dBx);
    elseif method == 2
        
    end