function x = qba(a, B, d, epssub)
[m, n] = size(d);
p = 100;
if a == 0
    [dm, imax] = max(d);
    et = zeros(m, 1);
    et(imax) = 1;
else
    F = 1 : n;
    ninf0 = 2 * n;
    while ninf0 >0
        Fold = F;
        [x1, F, iter, ninf, ~] = BBP(A, - d / a, F, p, ninf0, epssub, 0, debug);
        if iter == p
            [x1, F, iter, ninf, ~] = BBP(A, - d / a, Fold, p, ninf0, epssub, 1, debug);
        end
        ninf0 = ninf;
    end
    a = d' * x1/