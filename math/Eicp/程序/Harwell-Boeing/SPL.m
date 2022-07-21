%% SPL
function [x, iter]=SPL(A, B, x0, maxIT, eps, epsbbp)
[~, n] = size(A);
x = x0;
iter =0;
typeB = 0;
if isdiag( B )
   if isequal(B, eye(n))
    typeB = 2;
   else
    typeB = 1;
    B = diag(B);
   end
end
p = min(n, 20);
F = 1 : n;

while 1
    iter = iter + 1;
    Ax = A * x;
    Bx = computexB(B, x, typeB);
    xAx = x' * Ax;
    xBx = x' * Bx;
    lamdab = xAx / xBx;
    yk = -lamdab * Bx;
    ninf0 = 2 * n;
    while ninf0 >0
        Fold = F;
        [x1, F, iterk, ninf] = BBP(A, yk, F, p, ninf0, epsbbp, 0);
        if iterk == p
        [x1, F, ~, ninf] = BBP(A, yk, Fold, p, ninf0, epsbbp, 1);
        end
        ninf0 = ninf;
    end
    d = x1 - x;
    x = x1 ; 
    if norm( d ) <= eps || iter > maxIT
        break
    end
end
end
function Bx = computexB(B, x, typeB)
   if typeB == 1
        Bx  = B .* x;
   elseif typeB == 2
       Bx = x;
   else
        Bx = B * x;
   end
end
