%[xk,resvec,arvec,face1vec,face2vec,tf]
function [xk, resvec, arvec, faceXvec, tf] = IPG(A,b,x0,tol,det,tou,maxit)
%x
t=clock;
loopcount = 0;
[rpk, normr, ~, ~, normKKT, faceX, ~] = kktResidual(A, b, x0 , []);
Ab = A'*b;
I = rpk>0;
g = 
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
    if loopcount > maxit
        break;
    end
    loopcount = loopcount +1;
    I = rpk >0;
    Ax = b(I) - rpk(I);
    ADA = A(I,:)' * Ax;
    d = x0./(ADA+det);
    p = -d .* g;
    alphaAll = -x0./p;
    alphak = min(alphaAll(alphaAll>0));
    talphak = tou * alphak;
    [alpha, ~] = wolfe(A, b, x0, p, talphak, 20);
    %[alpha, knots, retcode] = arraySpiece(A,b,x0,p);
    x0 = x0 + talphak*p;
    [rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
    % record the value of objection function
    resvec(loopcount) = normr;
    % record the value of the gradient function
    arvec(loopcount) = normKKT;
    faceXvec(loopcount) = face1;
    
end
xk=x0;
tf = etime(clock,t);
end

