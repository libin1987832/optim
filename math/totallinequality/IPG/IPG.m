%[xk,resvec,arvec,face1vec,face2vec,tf]
function [xk,resvec,arvec,face1vec,face2vec,tf] = IPG(A,b,x0,tol,det,tou,maxit)
%x
t=clock;
loopcount = 1;
[rpk, normr, xmin, g, normKKT, face1, face2] = kktResidual(A, b, x0 , [], 1);
% the residual vector
resvec = zeros(1,maxit);
% the normal gradient
arvec = zeros(1,maxit);
% subspace minization
itersm = zeros(1,maxit);
% face b-Ax
face1vec = zeros(1,maxit);
% face x
face2vec = zeros(1,maxit);
resvec(1) = xmin;
arvec(1) = normKKT;
face1vec(1) =face1;
face2vec(1) =face2;

%while norm(x0.*g,inf)>e
while norm(x0.*g,inf) > tol || min(g)< -tol
    if loopcount > maxit
        break;
    end
    loopcount = loopcount +1;
    D = max(diag(b-A*x0),0);
    D(D>0) = 1;
    d = x0./(A'*D*A*x0+det);
    p = -d .* g;
    alphaAll = -x0./p;
    alphak = min(alphaAll(alphaAll>0));
    talphak = tou * alphak;
    %a=fsearcha(A,b,x0,p); %ʹ��ԭ������������
    %a=fwolfepowersearcha(A,b,x0,p);%ʹ��wolfepower��������������������
    %    [alpha,x0,fx0,g]=wolfe(A, b, x0, p, 0.9);
    [alpha, knots, retcode] = arraySpiece(A,b,x0,p);
    if alpha < talphak
        x0 = x0 + alpha*p;
    else
        x0 = x0 + talphak*p;
    end
    [rpk, normr, xmin, g, normKKT, face1, face2] = kktResidual(A, b, x0 , [], 1);
    % record the value of objection function
    resvec(loopcount) = normr;
    % record the value of the gradient function
    arvec(loopcount) = normKKT;
    face1vec(loopcount) = face1;
    face2vec(loopcount) = face2;
end
xk=x0;
tf = etime(clock,t);
end

