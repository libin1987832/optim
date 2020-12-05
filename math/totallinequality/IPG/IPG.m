%[xk,resvec,arvec,face1vec,face2vec,tf]
function [xk, resvec, arvec, faceXvec, tf] = IPG(A,b,x0,tol,det,tou,maxit,dtype)
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
            
    %cos1 = -g'*p / (norm(g)*norm(p));
    switch dtype
        case 'NT'
            p = -(A(I,:)'*A(I,:)+det)\g;
            %  cos2 =  np'*p / (norm(np)*norm(p));
        case 'ST'
            p = -g;
        otherwise
            Ax = b - rpk;
            ADA = A(I,:)' * Ax(I);
            d = x0./(ADA + det);
            if sum(abs(d<0))>0
            assert(sum(abs(d<0)) > 0);
            end
            p = - d .* g;
            assert( g'*p < 0);
            
    end
    x0(x0<1e-15) = 0;
    p(x0<1e-15) = 0;
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
%% test example
% A = [1 3;2 4;-5 -6]; b = [5;6;-3]; x0 = [1;1];
% -g = [39;45]; ADA = [26 33;33 45]; ADAx = [59;78];
% p = [39/59;45/78] pnt = [30/9; -13/9] 
% alpha = +inf alpha = 9/13  alpha =1 f = 8+39*5/59+180/78=13.6128^2

