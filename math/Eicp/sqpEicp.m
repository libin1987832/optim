function [x, iter, error] = sqpEicp(A, B, M, x0, sigma0, dsigma, eps, bisectEps, maxIT, debug)
error = 0;
x = x0;
[m, n] = size(M);
if debug
    error = zeros(1, maxIT);
end
iter = 0;
method = 1;
typeM = 0;
typeB = 0;
MA = zeros(n,1);
if isdiag(M)
    method = 2;
    typeM = 2;
    M = diag(M);
    MA = M - diag(A);
end
if isequal( A, M )
    typeM = 1;
end
if isequal(M, eye(n))
    typeM = 4;
end
if isequal(B, eye(n))
    typeB = 2;
end
if typeM==0 && method == 1 && isdiag( M - A )
    typeM = 3;
    MA = diag( M - A);
end
if isdiag( B )
    typeB = 1;
    B = diag(B);
end
sigma = 0;
while 1
    iter = iter + 1;
    xB = computexB(B, x, typeB)';
    Ax = A * x;
    tao = 1 - 0.5 * xB * x;  
    % compute the dirction
    [d, lambda, dAx, Md] = searchdir(M, Ax, xB, tao, x, method, bisectEps);
    if norm( d ) < eps || iter > maxIT
        break;
    end
    Ad = computeAd(A, d, Md, MA, typeM);
    dAd = d' * Ad;
    Bd = computexB(B, d, typeB);
    c = 0.5 * d' * Bd;
    e=abs(lambda)+dsigma;
    if e > sigma
        sigma = e;
    end
    z= - ( dAx + tao * sigma ) / ( dAd + 2 * c * sigma );
    u= - ( tao - sqrt( tao^2 + 4 * c * tao ) ) / ( 2 * c );
    % compute steplength
    alpha = searchstep(tao, z, u);
    x = x + alpha * d;
    if debug
        error(1,iter) = 0.5 * x' * Ax;
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
function Ad = computeAd(A, d, Md, MA, typeM)
    if typeM == 1
        Ad = Md;
    elseif typeM == 3
        Ad = Md - MA .* d;
    else
        Ad = A * d;
    end
end