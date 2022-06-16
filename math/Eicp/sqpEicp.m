function [x, error] = sqpEicp(A, B, M, x0, sigma0, eps, bisectEps, maxIT, debug)
error = 0;
x = x0;
if debug
    error = zeros(1, maxIT);
end
k = 0;
method = 1;
if norm(diag(D) - M) == 0
    method = 2;
end
maxe = sigma0;
while 1
    xB = x' * B;   
    Ax = A * x;
    tao = 1 - 0.5 * xB * x;  
    % compute the dirction
    d = searchdir(Ax, xB, tao, x, method, bisectEps);
    if norm( d ) < eps || k > maxIT
        break;
    end
    Ad = A * d;
    dAd = d' * Ad;
    Bd = B * d;
    c = 0.5 * d' * Bd;
    e=abs(lambda)+0.01;
    if e > maxe
        sigma = e;
    end
    z= - ( Ax' * d + tao * sigma ) / ( dAd + 2 * c * sigma );
    u= ( - tao + sqrt( tao^2 + 4 * c * tao ) ) / ( 2 * c );
    % compute steplength
    alpha = searchstep(tao, z, u);
    x = x + alpha * d;
    k = k + 1;
    if debug
        error(1,k) = 0.5 * x' * A * x;
    end
end