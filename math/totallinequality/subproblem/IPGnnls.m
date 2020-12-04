function [x0, res] = IPGnnls(A, b, x0, rpk, tol, tou, maxit)
Ax = b - rpk;
ATb = A' * b;
ATAx = A' * Ax;
Ar = ATAx - ATb;
loopcount = 1;
res = zeros(maxit,1);
while max(abs(x0.* Ar)) > tol || min( Ar ) > -tol
    d = x0./ATAx;
    p = - d .* Ar;
    alphaAll = - x0./p;
    alphak = min(alphaAll(alphaAll>0));
    if isempty(alphak)
        alphak = 1;
    end
    talphak = tou * alphak;
    Ap = A * p;
    Apn = Ap' * Ap;
    alphaStar = - (p' * Ar) / Apn;
    alpha = min(talphak, alphaStar);
    x0 = x0 + alpha * p;
    Ax = A * x0;
    ATAx = A' * Ax;
    Ar = ATAx - ATb;
    loopcount = loopcount +1;
    r = b - Ax;
    res(loopcount) = 0.5 * (r' * r); 
    if loopcount > maxit
        break;
    end
end

%     xa = [0:0.0001:talphak];
%     ya = arrayfun(@(alpha) func(A,b,x0,p,alpha), xa);
%     ya2 = arrayfun(@(alpha) normr + alpha * 0.25 * g' * p , xa);
%     ya3 = arrayfun(@(alpha) func(A,b,x0,p,alpha), allalpha);
%     plot(xa,ya,'+',xa,ya2,'o',allalpha,ya3,'x');
%     hold on

%[alpha, knots, retcode] = arraySpiece(A,b,x0,p);
%        x0 = x0 + alpha*p;
