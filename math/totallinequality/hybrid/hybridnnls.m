function [xk, resvec, arvec, face1vec, face2vec, tf] = hybridnnls(A, b, x0, tol, nf, maxit, options, type)
t=clock;
% stop criterion
display = false;
[m,n] = size(A);
[rpk, normr, ~, g, normKKT, face1, face2] = kktResidual(A, b, x0,[],1);
iter = 1;
maxit = 2 * maxit;
% the residual vector
resvec = zeros(1,maxit);
% the normal gradient
arvec = zeros(1,maxit);
% face b-Ax
face1vec = zeros(1,maxit);
% face x
face2vec = zeros(1,maxit);
resvec(1) = normr;
arvec(1) = normKKT;
face1vec(1) =face1;
face2vec(1) =face2;

% flag 0-4 return lsqr flag
flag = 5;
while norm( x0 .* g, inf) > tol || min( g )< -tol
    iter = iter + 2;
    [xkA, rpk] = simple(A, b, x0, n, rpk, nf, 100*tol, options, type);
    [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, xkA(:, end), [], 1);
    resvec(iter) = normr;
    % record the value of the gradient function
    arvec(iter) = normKKT;
    face1vec(iter) = face1;
    face2vec(iter) = face2;
    if display
        %pg = g'*p;
        [minx,loc]=min(x0);
%         fprintf('simple(%d end): norm(%g),normKKT(%g),minx(%g,%d),gp(%g),xa(%d),ba(%d)\n',...
%             (iter-1)/2,normr,normKKT,minx,loc,pg,face1,face2);
%        fprintf('simple(%d end): norm(%g),normKKT(%g),minx(%g),xa(%d),ba(%d)\n',(iter-1)/2,normr,normKKT,minx,face1,face2);
         fprintf('simple(%d end): norm(%g),normKKT(%g),minx(%g,%d),xa(%d),ba(%d)\n',...
             (iter-1)/2,normr,normKKT,minx,loc,face1,face2);
    end
      %  isSub = true;
    isSub = strategy(A,b,x0,[],'PHA',iter,nf,rpk,xkA);
    if isSub
        [x0, rpk] = newton(A, b, n, rpk, xkA(:, end));
        [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, x0 , rpk, 1);
        % record the value of objection function
        resvec(iter +1) = normr;
        % record the value of the gradient function
        arvec(iter +1) = normKKT;
        face1vec(iter +1) = face1;
        face2vec(iter +1) = face2;
        if display
            pg = g'*(x0-xkA(:, end));
            [minx,loc]=min(x0);
            fprintf('newton(%d end): norm(%g),normKKT(%g),minx(%g,%d),gp(%g),xa(%d),ba(%d)\n',...
                (iter-1)/2,normr,normKKT,minx,loc,pg,face1,face2);
        end
    else
        x0 = xkA(:, end);
    end
    
    %    if iter > maxit || flag == 0
    if iter > maxit
        break;
    end
end
xk = x0;
% resvec = resvec(resvec>0);
% % the normal gradient
% arvec = arvec(arvec>0);
% % face b-Ax
% face1vec = face1vec(face1vec>0);
% % face x
% face2vec = face2vec(face2vec>0);
tf = etime(clock,t);


%    [alpha, minf, knot] = arraySpiecewise(A,b,x0,p);
%
% %         knoty = arrayfun(@(alpha) funmin(A,b,x0,u,alpha), knot);
%
%         xa = [0.019:0.00001:0.021];
% ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
%  pxy={};
% % pxy(1).X = knot;
% % pxy(1).Y = knoty;
% pxy(1).X = xa;
% pxy(1).Y = ya;
% figure
% hold on
% p1 = arrayfun(@(a) plot(a.X,a.Y),pxy);
% p1(1).Marker = 'o';
% %p1(2).Marker = '+';
% hold off
% r = b - A * x0;
% r(r<0)=0;
% Ar = A'*r;
% [alpha,Ar'*p]
