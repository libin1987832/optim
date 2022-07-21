function x = FLQP(x0, B, A, n, strategy, maxIt, maxItsub, eps, epsbbp)
Bx0 = B * x0;
Ax0 = A * x0;
xAx = x0' * Ax0;
xBx = x0' * Bx0;
a = xAx / ( 2 * Bx0' * x0 - xBx);
a0 = a + 10;
iter = 0;

 while abs(a - a0) > eps && iter < maxIt % && abs(a) > eps
a0 = a;
iter = iter + 1;

if strategy == 1
    [x, ~, ~] = BBP2(Bx0, A, n, a, 0.5 * xBx, maxItsub, epsbbp);
else
    opts = optimset('Display','off');
    x = quadprog(A, -a * Bx0, -Bx0', -0.5*xBx , ones(1, n), 1, zeros(n, 1), [ ], [], opts);
end
 Ax0 = A * x;
 xAx = x' * Ax0;
 a = xAx / ( 2 * Bx0' * x - xBx);
 end


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


