%% SPL（别人的算法）
function [x, eta, iters, nitBBs, error]=SPL(A, B, x0, maxIT, eps, epssub, debug)
[m, n] = size(A);
x = x0;
iters =0;
typeB = 0;
if size( B, 2 ) == 1
    typeB = 1;
end
error = 0;
if debug
    error = zeros(maxIT, 1);
end
p = min(n, 20);
nitBBs = 0;
N = 1 : n;
while 1
    Ax = A * x;
    Bx = computexBx(B, x, typeB);
    xAx = x' * Ax;
    xBx = x' * Bx;
    lamdab = xAx / xBx;
    yk = -lamdab * Bx;
    if iters == 0
        F = N ;
        d = eps + 1;
    end
    if debug
        error(iters+1) =  xAx / xBx;
    end
    w = Ax - lamdab * Bx;
    if norm( d ) <= eps || iters >= maxIT || min(w) >= -eps
        break
    end
   T = setdiff( N , F );
    ninf0 = sum( x < 0 ) + sum( w < 0  );
    while ninf0 >0
        Fold = F;
        [xeta1, F, nitBB, ninf, ~] = BBP(A, yk, F, p, ninf0, epssub, 0, debug);
        if nitBB == p
            [xeta1, F, nitBB, ninf, ~] = BBP(A, yk, Fold, p, ninf0, epssub, 1, debug);
        end
        nitBBs = nitBBs + nitBB;
        ninf0 = ninf;
    end
    x1 = xeta1(1:n);
    eta = xeta1(n+1);
    d = x1 - x;
    x = x1 ; %迭代更新
    iters = iters + 1;
end
nitBBs = nitBBs / iters;
end
function Bx = computexBx(B, x, typeB)
if typeB
    Bx  = B .* x;
else
    Bx = B * x;
end
end