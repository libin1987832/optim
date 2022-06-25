function [x, iter, error] = sqpEicp(A, B, M, x0, sigma0, dsigma, eps, bisectEps, maxIT, debug)
error = 0;
x = x0;
if debug    
    error = zeros(1, maxIT);
end
iter = 0;
method = 1;
typeM = 0;
typeB = 0;
if isdiag(M)
    method = 2;
    typeM = 2;
    M = diag(M); 
end
if isequal( A, M )
    typeM = 1;
end 
if typeM==0 && method == 1 && isdiag( M - A ) % 不明白矩阵M不就是分两种吗？
    typeM = 3;
    
end
if isdiag( B )
    typeB = 1;
end
sigma = sigma0;
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
    Ad = computeAd(A, d, Md,  typeM);
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
   if typeB
        Bx  =diag(B) .* x;
    else
        Bx = B * x;
   end
end
function Ad = computeAd(A, d, Md, typeM) 
    if typeM == 1
        Ad = Md;
    elseif typeM == 3   %不明白这是为什么   
        MA = diag( M - A);
        Ad = Md - MA .* d;
    else
        Ad = A * d;
    end
end