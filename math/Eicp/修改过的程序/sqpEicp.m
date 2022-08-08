function [x, iter] = sqpEicp(A, B, M, x0, sigma0, dsigma, eps, bisectEps, maxIT, typeM, method)

x = x0;
iter = 0;
sigma = sigma0;
if method == 2
  M = diag(M);
end    
while 1
    iter = iter + 1;
    xB = computexB(B, x)';
    Ax = A * x;
    tao = 1 - 0.5 * xB * x;  
    % compute the dirction
    
    [d, lambda, dAx, Md] = searchdir(M, Ax, xB, tao, x, method, bisectEps);
    if norm( d ) <= eps || iter > maxIT
        break;
    end
    Ad = computeAd(A, d, M, Md,  typeM);
    dAd = d' * Ad;
    Bd = computexB(B, d);
    c = d' * Bd;
    e=abs(lambda)+dsigma;
    
    if e > sigma
        sigma = e;
    end

    z= - ( dAx + tao * sigma ) / ( dAd + c * sigma );
    u= - ( tao - sqrt( tao^2 + 2 * c * tao ) ) / c ;
    % compute steplength
    alpha = searchstep(tao, z, u);
    x = x + alpha * d;
end
end

function Bx = computexB(B, x)
   if isdiag( B )
        Bx  =diag(B) .* x;
    else
        Bx = B * x;
   end
end

function Ad = computeAd(A, d, M, Md, typeM) 

    if typeM == 1
        Ad = Md; 
    elseif typeM == 2   
        MA = diag( M - A);
        Ad = Md - MA .* d;
    else
        Ad = A * d;
    end
end