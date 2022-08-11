%% SPL�����˵��㷨��
function [x, crit, iters, nitBBs, error]=SPL(A, B, x0, maxIT, eps1, eps2, epssub, strategy,debug)
[m, n] = size(A);
x = x0;
iters =0;
typeB = 1;
if size( B, 2 ) == 1
    typeB = 2;
end
if size( B, 2 ) == n && (isequal(speye(n), B) || isequal(ones(n,n), B))
    typeB = 3;
end
error = 0;
if debug
    error = zeros(maxIT, 1);
end
p = min(n, 20);
nitBBs = 0;
N = 1 : n;
yk = 0;
eta3 = 0;
while 1
    Ax = A * x;
    Bx = computexBx(B, x, typeB);
    xAx = x' * Ax;
    xBx = x' * Bx;
    lamdab = xAx / xBx;
    %eta = Ax + yk - eta3*ones(n,1);
    eta = Ax - lamdab * Bx;
    if min(eta) >= -eps2 
        crit = 2;
        break;
    end

    if iters == 0
        F = N ;
    end
    if debug
        error(iters+1) =  xBx / xAx;
    end
    yk = -lamdab * Bx;

   if strategy == 0
        [x1, ~, F, nitBB, ~] = BBP3(A, yk, F, 20, epssub, debug);
        nitBBs = nitBBs + nitBB;
   elseif strategy == 1
        ninf0 = sum( x < 0 ) + sum( eta < 0  );
        while ninf0 >0
            Fold = F;
            [x1, F, nitBB, ninf, ~] = BBP(A, yk, F, p, ninf0, epssub, 0, debug);
            if nitBB == p
                [x1, F, nitBB, ninf, ~] = BBP(A, yk, Fold, p, ninf0, epssub, 1, debug);
            end
            nitBBs = nitBBs + nitBB;
            ninf0 = ninf;
        end
   else
       opts = optimoptions('quadprog','Display','off');
       [x1,FVAL,EXITFLAG,OUTPUT,LAMBDA] =quadprog(A, yk, [], [], ones(1,n), 1, zeros(n,1), Inf*ones(n,1),x, opts);
        nitBBs = nitBBs + 1;
   end
    d = x1 - x;
    x = x1 ; %��������
    iters = iters + 1;

    if norm( d ) <= eps1
        crit = 1;
        break;
    end
    
    if iters >= maxIT
        crit = 3;
        break;
    end
end
nitBBs = nitBBs / iters;
end
function Bx = computexBx(B, x, typeB)
if typeB == 2
    Bx  = B .* x;
elseif typeB == 3
    Bx = x;
else
    Bx = B * x;
end
end