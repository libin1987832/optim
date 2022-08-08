function [x, iter] = qba(x0, B, d, eps, epssub)
[m, n] = size(d);
p = 100;
Ax0 = A * x0;
x0Ax0 = x0 * Ax0;
a = d' * x0 / x0Ax0 ;
if a == 0
    [~, imax] = max(d);
    x = zeros(m, 1);
    x(imax) = 1;
end
while 1
    F = 1 : n;
    ninf0 = 2 * n;
    while ninf0 >0
        Fold = F;
        [x1, F, iter, ninf, ~] = BBP(B, - d / a, F, p, ninf0, epssub, 0, debug);
        if iter == p
            [x1, F, iter, ninf, ~] = BBP(B, - d / a, Fold, p, ninf0, epssub, 1, debug);
        end
        ninf0 = ninf;
    end
    Bx = B * x1;
    xBx = x * Bx;
    a1 = (d' * x1) / xBx;
    if abs( a - a1 ) >= eps
        x = x1;
    else
        a = a1;
    end
end