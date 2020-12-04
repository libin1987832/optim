% 网上下载的一个wolf条件实现算法
function [alpha, retcode] = wolfe(A,b,xk, dk, range, maxit)
% alphaMin=t*min(-1.*xk./dk);
rho = 0.25; sigma = 0.75;
[~, normr, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk , [], []);
qpc1 = rho * Ar' * dk;
qpc2 = sigma * Ar' * dk;
alpha = range;
loopcount = 0;
while func(A , b, xk, dk ,alpha) > normr +  alpha * qpc1   
    loopcount = loopcount + 1;    
    alpha=alpha/2;
    if loopcount > maxit
        alpha = range;
        retcode = [2,loopcount];
        return;
    end
end
sumloop = loopcount; 
loopcount = 0;
beta = 2 * alpha;
[~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + alpha * dk , [], []);
while Ar' * dk < qpc2 
    loopcount = loopcount + 1;    
    alpha=alpha/2;
    [~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + alpha * dk , [], []);
    if loopcount > maxit
        alpha = range;
        retcode = [3,sumloop];
        return;
    end
end
middle = 0.5 * (alpha + beta);
[~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + middle * dk , [], []);
Ardk = Ar' * dk; 
sumloop = sumloop + loopcount; 
loopcount = 0;
while Ardk < qpc2 
    loopcount = loopcount + 1;    
    if Ardk > 0
    beta = middle;
    else
    alpha = middle;    
    end
    middle = 0.5 * (alpha + beta);
    [~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + middle * dk , [], []);
    Ardk = Ar' * dk; 
    if loopcount > maxit
        alpha = range;
        retcode = [4,sumloop];
        return;
    end
end
retcode = [1,sumloop + loopcount];
end 
function fvalue = func(A,b,x0,p,alpha)
r = b - A * (x0 + alpha * p);
r( r < 0 ) = 0;
fvalue = 0.5 * (r' * r);
end
