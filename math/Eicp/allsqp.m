function [x, iter, error] = allsqp(A, B, M, x0, sigma0, dsigma, eps, bisectEps, maxIT, debug)
n = size(A, 2);
np = 2^n - 1;
x = zeros(n, np);
iter = zeros(1,np);
error =zeros(1,np);
charA = char(ones(n,1)*double('1'))';
for iter_i = 1 : np
    F = dec2bin(iter_i, n) == charA;
    AF = A(F, F);
    BF = B(F, F);
    MF = M(F, F);
    xF = x0(F);
    O=inv(BF)*AF;
    R=max(abs(eig(O))); %A���װ뾶  
    sigma0 = R + 0.01;
    [xmin,  itermin, errormin] = sqpEicp(AF, BF, MF, xF, sigma0, dsigma, eps, bisectEps, maxIT, debug);
    [xmax,  itermax, errormax] = sqpEicp(BF, AF, MF, xF, sigma0, dsigma, eps, bisectEps, maxIT, debug);
    if sum(xmin > 0) == sum(F)
        x1 = xmin;
        iter1 = itermin;
        error1 = errormin;
    else
        x1 = xmax;
        iter1 = itermax;
        error1 = errormax;
    end
    x(F, iter_i) = x1;
    iter(1, iter_i) = iter1;
    error(1, iter_i) = error1;
end