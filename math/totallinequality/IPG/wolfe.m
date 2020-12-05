% 网上下载的一个wolf条件实现算法
function [alpha, allalpha, retcode] = wolfe(A,b,xk, dk, range, maxit)
% alphaMin=t*min(-1.*xk./dk);
rho = 0.001; sigma = 0.8;
%[~, normr, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk , [], []);
[normr,Ar,~,~] = dffunc(A,b,xk);
qpc1 = rho * Ar' * dk;
qpc2 = sigma * Ar' * dk;
alpha = range;
loopcount = 1;
allalpha = zeros(3 * maxit ,1);
allalpha(1) = alpha; 
while func(A , b, xk, dk ,alpha) > normr +  alpha * qpc1   
    loopcount = loopcount + 1;    
    alpha=alpha/2;
    allalpha(loopcount) = alpha; 
    if loopcount > maxit
        alpha = range;
        retcode = [2,loopcount];
        return;
    end
end
sumloop = loopcount; 
loopcount = 0;
beta = 2 * alpha;
% [~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + alpha * dk , [], []);
while Ar' * dk < qpc2 
    loopcount = loopcount + 1;    
    alpha=alpha/2;
    allalpha(sumloop + loopcount) = alpha; 
%     [~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + alpha * dk , [], []);
    Ar = dfunc(A,b,xk,dk,alpha);    
    if loopcount > maxit
        alpha = range;
        retcode = [3,sumloop];
        return;
    end
end
middle = 0.5 * (alpha + beta);
% [~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + middle * dk , [], []);
Ar = dfunc(A,b,xk,dk,middle);    
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
    allalpha(sumloop + loopcount) = alpha; 
    %[~, ~, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk + middle * dk , [], []);
    Ar = dfunc(A,b,xk,dk,middle);    
    Ardk = Ar' * dk; 
    if loopcount > maxit
        alpha = range;
        retcode = [4,sumloop];
        return;
    end
end
retcode = [1,sumloop + loopcount];

%     xa = [0:range/50:range]; 
%     ya = arrayfun(@(alpha) func(A,b,xk,dk,alpha), xa);   
%     ya2 = arrayfun(@(alpha) normr + alpha * qpc1 , xa);
%     ya4 = arrayfun(@(alpha) dfunc(A,b,xk,dk,alpha) >= qpc2 , xa);
%     ya3 = arrayfun(@(alpha) func(A,b,xk,dk,alpha), allalpha);
%     plot(alpha,func(A,b,xk,dk,alpha),'v',...
%         xa,ya,'+',xa,ya2,'o',allalpha,ya3,'x',...
%         xa(ya4),ones(size(xa(ya4)))*ya(1),'-');
%     hold on

end 
function fvalue = func(A,b,x0,p,alpha)
r = b - A * (x0 + alpha * p);
r( r < 0 ) = 0;
fvalue = 0.5 * (r' * r);
end
function fvalue = dfunc(A,b,x0,p,alpha)
Ap = A * p;
r = b - A * (x0 + alpha * p) ;
r( r < 0 ) = 0;
fvalue =  -Ap' * r;
end
function [fvalue,gvalue,normKKT,Ax] = dffunc(A,b,x0)
Ax = A * x0;
r = b - Ax;
r( r < 0 ) = 0;
gvalue =  -A' * r;
normKKT = norm( x0 .* gvalue, inf);
fvalue = 0.5 * (r' * r);
end