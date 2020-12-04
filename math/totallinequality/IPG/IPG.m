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
%      if normr <40.9262
%          dd=1;
%      end
    if loopcount > maxit
        break;
    end
    loopcount = loopcount +1;
    I = rpk > 0;
    Ax = b - rpk;
    ADA = A(I,:)' * Ax(I);
    d = x0./(ADA+det);
    assert(sum(abs(d<0))==0);
    p = - d .* g;
    assert( g'*p < 0);
    cos1 = g'*p / (norm(g)*norm(p));
    np = -(A(I,:)'*A(I,:))\g;  
    cos2 =  np'*p / (norm(np)*norm(p))
    
    alphaAll = - x0./p;
    alphak = min(alphaAll(alphaAll>0));
    if isempty(alphak)
        alphak = 10;
    end
    talphak = tou * alphak;

    [alpha, ~, ~] = wolfe(A, b, x0, p, talphak, 20);
    x0 = x0 + alpha * p;
    [rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
    % record the value of objection function
    resvec(loopcount) = normr;
    % record the value of the gradient function
    arvec(loopcount) = normKKT;
    faceXvec(loopcount) = faceX;
end
xk=x0;
tf = etime(clock,t);
end

