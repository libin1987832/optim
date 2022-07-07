%% SPL（别人的算法）
function [x, iters, error]=SPL(A, B, x0, maxIT, eps, epssub, debug)
[m, n] = size(A);
x = x0;
iters =0;
typeB = 0;
if size( B, 2 ) == 1
    typeB = 1;
end
error = 0;
if debug
    error = zeros(maxIts, 1);
end
p = min(n, 20);
while 1
    iters = iters + 1;
    Ax = A * x;
    Bx = computexBx(B, x, typeB);
    xAx = x' * Ax;
    xBx = x' * Bx;
    lamdab = xAx / xBx;
    yk = -lamdab * Bx;
    if iters == 1
        F = 1 : n;
    end
    ninf0 = 2 * n;
    while ninf0 >0
        Fold = F;
        [x1, F, iter, ninf, ~] = BBP(A, yk, F, p, ninf0, epssub, 0, debug);
        if iter == p
            [x1, F, iter, ninf, ~] = BBP(A, yk, Fold, p, ninf0, epssub, 1, debug);
        end
        ninf0 = ninf;
    end
    d = x1 - x;
    x = x1 ; %迭代更新
    if debug
        error(iters) = - xAx / xBx;
    end
    if norm( d ) <= eps || iters > maxIT
        break
    end
end
end
function Bx = computexBx(B, x, typeB)
if typeB
    Bx  = B .* x;
else
    Bx = B * x;
end
end