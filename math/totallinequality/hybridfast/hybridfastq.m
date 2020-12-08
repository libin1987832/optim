function [xk, resvec, arvec, face1vec, face2vec, tf] = hybridfast(A, b, x0, tol, nf, maxit)
t=clock;
% stop criterion
display = true;
if display 
    maxit = 2;
end
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
AT = A';
% flag 0-4 return lsqr flag
flag = 5;
p = -g;
lsqrTol= 1e-5;
maxIter =10;
maxKnot = 1e5;
testalphax = zeros(20);
testalphab = zeros(20);
maxits = 10;
eIter = 3;
while norm( x0 .* g, inf) > tol || min( g )< -tol
    iter = iter + 2;
    for i = 2:nf+1
        alphaAll = - x0./p;
        alphak = man(alphaAll(alphaAll>0));
        steplength = alphak;
        rho = 0.1;
        qpc1 = - rho * p' * dk;
        alpha = steplength;
        loopcount = 1;
        allalpha = zeros(3 * maxit ,1);
        allalpha(1) = alpha;
        while func(A, b, Ax, Adk ,alpha) > normr +  alpha * qpc1
            loopcount = loopcount + 1;
            alpha=alpha/2;
            allalpha(loopcount) = alpha;
            if loopcount > maxit
                alpha = range;
                retcode = [2,loopcount];
                return;
            end
        end
        retcode = [1,loopcount];
        return
        end
        if display
            [alpha, minf, knot,retcode] = arraySpiece(A,b,x0,p,tol,maxits);
            [~, normr, ~, g, normKKT, face1, faceN] = kktResidual(A, b, x0, [], 1);
            xt = x0 + steplength * p; xt(xt<0) =0;
            [~, normrN, ~, g0, normKKT, faceN1, faceN2] = kktResidual(A, b, xt, [], 1);
            pg = g0' * -g; 
            fprintf('simple(steepdown,%d): normr(%g,%g),pg(%g),activeX(%d,%d),activeB(%d,%d) alpha(%g,%g)\n'...
                ,i,normr,normrN,pg,face1,faceN1,face2,faceN2,steplength,alpha);
            knoty = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), [testalphax,testalphab]);
        end
        x0 = x0 + steplength * p;
        x0(x0<0) =0;
        if i < nf+1
            r = b - A * x0;
            I = r > 1e-10;
            p = AT(: , I) * r(I);
        end
    end
    
    if display
        [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, x0, [], 1);
        fprintf('simple(%d end): norm(%g),normKKT(%g),minx(%g),xa(%d),ba(%d)\n',...
            (iter-1)/2,normr,normKKT,minx,face1,face2);
        resvec(iter) = normr;
        % record the value of the gradient function
        arvec(iter) = normKKT;
        face1vec(iter) = face1;
        face2vec(iter) = face2;
    end
    
%     rkn = r - eIter * A * p;
    rkn = r - eIter * ApIu;
    Irkn = rkn > tol;
    ssign=sum(~xor(I , Irkn));
    
    if ssign==m
        RR = x0 > 1e-5;
        bI = r(I);
        AI = A(I,RR);
        [u,flag,relres,iterlsqr,resvec,lsvec,out] = lsqrmx( AI ,bI,lsqrTol,maxIter,[],[],zeros(size(AI,2),1),A(:,RR),b,x0(RR),I,3);
        
        if norm(u) > 1e-10
            x0(RR) = x0(RR) + u;
        end
        if display
            [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, x0 , rpk, 1);
            % record the value of objection function
            resvec(iter +1) = normr;
            % record the value of the gradient function
            arvec(iter +1) = normKKT;
            face1vec(iter +1) = face1;
            face2vec(iter +1) = face2;
            pg = g(RR)' * u;
            [minx,loc]=min(x0);
            x0 + u;
            fprintf('newton(%d end): norm(%g),normKKT(%g),minx(%g,%d),gp(%g,%g),xa(%d,%d)\n',...
                (iter-1)/2,normr,normKKT,minx,loc,pg,norm(u),face1,face2);
        end
    end
    
    
    r = b - A * x0;
    I = r > 1e-10;
    p = AT(: , I) * r(I);
    %    if iter > maxit || flag == 0
    if iter > maxit
        break;
    end
end
xk = x0;
tf = etime(clock,t);
end
function f = funmin(A,b,x0,u,alpha)
xn = x0 + alpha*u;
xn( xn < 0) = 0;
f = b - A * xn;
f(f<0) = 0;
f = sqrt(f'*f);
end

function fvalue = func(A,b,Ax,Ap,alpha)
r = b - Ax - alpha * Ap;
r( r < 0 ) = 0;
fvalue = 0.5 * (r' * r);
end

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

