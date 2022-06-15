function [x, error] = sqpEicp(A, B, M, x0, method, eps, maxIT, debug)
error = 0;
x = x0;
if debug
    error = zeros(1, maxIT);
end
k = 0;
D = diag(M);
while 1
    xB = x' * B;   
    Ax = A * x;
    tao = 1 - 0.5 * xB * x;  
    d = searchdir(Ax, xB, tao, D, x, method);
    if norm( d ) < eps || k > maxIT
        break;
    end
    Ad = A * d;
    Bd = B * d;
    c = 0.5 * d' * Bd;
    
    sigma = 
    alpha = searchstep(tao, z, u);
    x = x + alpha * d;
    k = k + 1;
    if debug
        error(1,k) = 0.5 * x' * A * x;
    end
end