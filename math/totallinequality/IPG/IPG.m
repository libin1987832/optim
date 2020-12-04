%[xk,resvec,arvec,face1vec,face2vec,tf]
function [xk, resvec, arvec, faceXvec, tf] = IPG(A,b,x0,tol,det,tou,maxit)
%x
t=clock;
loopcount = 1;
[rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
% the residual vector
resvec = zeros(1,maxit);
% the normal gradient
arvec = zeros(1,maxit);
% face b-Ax
faceXvec = zeros(1,maxit);
resvec(1) = normr;
arvec(1) = normKKT;
faceXvec(1) =faceX;

while norm( x0 .* g, inf) > tol || min( g )< -tol
     if normr <40.9262
         dd=1;
     end
    if loopcount > maxit
        break;
    end
    loopcount = loopcount +1;
    %I = rpk >0;
    %Ax = b(I) - rpk(I);
    Ax = A * x0;
    I = Ax <= b;
    ADA = A(I,:)' * Ax(I);
    d = x0./(ADA+det);
    p = - d .* g;
    alphaAll = - x0./p;
    alphak = min(alphaAll(alphaAll>0));
    if isempty(alphak)
        alphak = 10;
    end
    talphak = tou * alphak;

    [alpha, allalpha, ~] = wolfe(A, b, x0, p, talphak, 20);
    
    xa = [0:talphak/50:talphak]; 
    ya = arrayfun(@(alpha) func(A,b,x0,p,alpha), xa);   
    ya2 = arrayfun(@(alpha) normr + alpha * 0.5 * g' * p , xa);
    ya3 = arrayfun(@(alpha) func(A,b,x0,p,alpha), allalpha);
    plot(xa,ya,'+',xa,ya2,'o',allalpha,ya3,'x');
    hold on
%     
    %[alpha, knots, retcode] = arraySpiece(A,b,x0,p);
    [rpk1, normr1, ~, g1, normKKT1, faceX1, ~] = kktResidual(A, b, x0 , [], 1);
    x0 = x0 + alpha*p;
    [rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
    if normr1 < normr
        norm1
    end
    %xa = [0.019:0.00001:0.021]; 
    %ya = arrayfun(@(alpha) func(A,b,x0,p,alpha), xa);

    % record the value of objection function
    resvec(loopcount) = normr;
    % record the value of the gradient function
    arvec(loopcount) = normKKT;
    faceXvec(loopcount) = faceX;
    
end
xk=x0;
tf = etime(clock,t);
end
function fvalue = func(A,b,x0,p,alpha)
r = b - A * (x0 + alpha * p);
r( r < 0 ) = 0;
fvalue = 0.5 * (r' * r);
end

