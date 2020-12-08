function [xk, resvec, arvec, face1vec, face2vec, tf] = hybridfast(A, b, x0, tol, nf, maxit)
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
AT = A';
% flag 0-4 return lsqr flag
flag = 5;
p = -g;
lsqrTol= 1e-5;
maxIter =10;
maxKnot = 1e5;
while norm( x0 .* g, inf) > tol || min( g )< -tol
    iter = iter + 2;
    for i = 2:nf+1
%         [~,x0, knot, retcode] = arraySpiece(A,b,x0,p,1e-5,30);
        nonzerou = p < 0;
        knot = -x0(nonzerou) ./ p(nonzerou);
        knot = sort(knot( knot > tol & knot < maxKnot));
        knot = [0 ; knot ; maxKnot];
        knot = unique(knot);
        left = 1;
        right = length(knot);
        loopcount = 1;
        while left+1  <= right && loopcount < maxits
            loopcount = loopcount + 1;
            knoti = floor(0.5*(left + right));
            x = x0 + 0.5*(knot(knoti) + knot(knoti+1))*p;
            Iu = x > 1e-10;
            [alpha,retcode] = spiecewise(A(:,~Iu),b,p(~Iu),x0(~Iu),knot( knoti ) , knot( knoti + 1 ));
            
            rIu = b - A(:,Iu) * x0(Iu);
            ApIu = A(:,Iu) * p(Iu);
            knotr = rIu ./ ApIu;
            knotr = sort( knotr( knotr >= knot( knoti ) & knotr <= knot( knoti + 1 )));
            knotr= unique([knot( knoti ); knotr ; knot( knoti + 1 )]);
for i = 1:length(knot)-1
    alpha = knot(i);
    beta = knot(i+1);
    % middle point give the active set
    rknot = r - 0.5 * (alpha + beta) * Ap;
    % Ap
    Apr = Ap(rknot > 0);
    % the derive Ap'(r-alpha*Ap)=0
    Ar = Apr' * r(rknot > 0);
    % alpha = (Ap*r)^T/(Ap'*Ap)
    steplength = Ar / ( Apr' * Apr );
    % the right point
    if steplength >= beta
        steplength = beta;
        retcode = [2,1];
    elseif steplength <= alpha
        steplength = alpha;
        retcode = [1,1];
        break;
    else
        retcode = [1,0];
        break;
    end
end
            
            
            
            if retcode(1) == 1 && retcode(2) == 0
                break;
            end
            if retcode(1) == 2 && retcode(2) == 0
                break;
            end
            if retcode(1) == 1 && retcode(2) ==1
                right = i;
            end
            if retcode(1) == 2 && retcode(2) ==1
                left = i;
            end
        end
        x0 = x0 + alpha * p;
        x0(x0<0) =0;  
        r = b - A * x0;
        I = r > 1e-10;
        p = AT(: , I) * r(I);
        if display
            [~, normrN, ~, g0, normKKT, faceN1, faceN2] = kktResidual(A, b, x0, [], 1);
            pg = g0' * -g; sumzeros = sum(xs > 1e-10);
            fprintf('simple(steepdown): normr(%g,),pg(%g),activeX(%d,%d),activeB(%d,%d)\n'...
                ,normr,normrN,pg,face1,faceN1,face2,faceN2);
        end
    end
    
    if display
        [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, xkA(:, end), [], 1);
        fprintf('simple(%d end): norm(%g),normKKT(%g),minx(%g,%d),xa(%d),ba(%d)\n',...
            (iter-1)/2,normr,normKKT,minx,loc,face1,face2);
        resvec(iter) = normr;
        % record the value of the gradient function
        arvec(iter) = normKKT;
        face1vec(iter) = face1;
        face2vec(iter) = face2;
    end
    
    rkn = r - eIter * A * p;
    Irkn = rkn > tol;
    ssign=sum(~xor(I , Irkn));
    
    if ssign==m
        RR = x0 > 1e-5;
        [u,flag,relres,iter,resvec,lsvec,out] = lsqrmx(A(I,RR),bI,lsqrTol,maxIter,[],[],zeros(size(AI,2),1),A(:,RR),b,x0(RR),AA,3);
        if display
            [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, x0 , rpk, 1);
            % record the value of objection function
            resvec(iter +1) = normr;
            % record the value of the gradient function
            arvec(iter +1) = normKKT;
            face1vec(iter +1) = face1;
            face2vec(iter +1) = face2;
            pg = g'*(x0-xkA(:, end));
            [minx,loc]=min(x0);
            fprintf('newton(%d end): norm(%g),normKKT(%g),minx(%g,%d),gp(%g),xa(%d),ba(%d)\n',...
                (iter-1)/2,normr,normKKT,minx,loc,pg,face1,face2);
        end
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
