function [x, error] = sqpEicp(A, B, M, x0, sigma0, eps, bisectEps, maxIT, debug)
error = 0;
x = x0;
if debug
    error = zeros(1, maxIT);
end
k = 0;
method = 1;
typeM = 0;
typeB = 0;
if isdiag(M)
    method = 2;
    M = diag(M);
end
if isequal( A, M )
    typeM = 1;
end 
if isdiag( M - A )
    typeM = 2;
    MA = diag( M - A );
end
if isdiag( B )
    typeB = 1;
    DB = diag(B);
end
maxe = sigma0;
while 1
    if typeB
        xB  = x .* B;
    else
        xB = x' * B;
    end
    Ax = A * x;
    tao = 1 - 0.5 * xB * x;  
    % compute the dirction
    [d, lambda, dAx, Md] = searchdir(M, Ax, xB, tao, x, method, bisectEps);
    if norm( d ) < eps || k > maxIT
        break;
    end
    if typeM == 1
        Ad = Md;
    elseif typeM == 2
        Ad = Md - MA .* d;
    else
        Ad = A * d;
    end
    dAd = d' * Ad;
    Bd = B * d;
    c = 0.5 * d' * Bd;
    e=abs(lambda)+0.01;
    if e > maxe
        sigma = e;
    end
    z= - ( dAx + tao * sigma ) / ( dAd + 2 * c * sigma );
    u= - ( tao - sqrt( tao^2 + 4 * c * tao ) ) / ( 2 * c );
    % compute steplength
    alpha = searchstep(tao, z, u);
    x = x + alpha * d;
    k = k + 1;
    if debug
        error(1,k) = 0.5 * x' * A * x;
    end
end