function [xk, resvec, arvec, face1vec, face2vec, tf] = hybridfast(A, b, x0, tol, nf, maxit)
t=clock;
% stop criterion
display = true;
if display 
    maxit = 100;
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
 xz = zeros(n,1);
 newton = 'ml';
 iterlsqr = 1;
while norm( x0 .* g, inf) > tol || min( g )< -tol
    iter = iter + 2;
    for i = 2:nf+1
        nonzerou = p < 0;
        knot = -x0(nonzerou) ./ p(nonzerou);
        knot = sort(knot( knot > tol & knot < maxKnot));
        knot = [0 ; knot ; maxKnot];
        knot = unique(knot);
        left = 1;
        right = length(knot);
        loopcountX = 0;
        while left+1  <= right && loopcountX < maxits
            loopcountX = loopcountX + 1;
            knoti = floor(0.5*(left + right));
            
            testalphax(loopcountX) = knot(knoti);
            
            x = x0 + 0.5*(knot(knoti) + knot(knoti+1))*p;
            Iu = x > 1e-10;
      
            rIu = b - A(:,Iu) * x0(Iu);
            ApIu = A(:,Iu) * p(Iu);
            knotr = rIu ./ ApIu;
            knotr = sort( knotr( knotr >= knot( knoti ) & knotr <= knot( knoti + 1 )));
            knotr= unique([knot( knoti ); knotr ; knot( knoti + 1 )]);
            leftb = 1;
            rightb = length(knotr);
            loopcountB = 0;
            while leftb + 1  <= rightb && loopcountB < maxits
                loopcountB = loopcountB + 1;
                knotri = floor(0.5*(leftb + rightb));
                alpha = knotr(knotri);
                beta = knotr(knotri + 1);
                
                testalphab(loopcountB) = alpha;
                
                % middle point give the active set
                rknot = rIu  - 0.5 * (knotr(knotri)+ knotr(knotri + 1)) * ApIu;
                % the derive Ap'(r-alpha*Ap)=0
                Apr = ApIu(rknot > 0);
                Ar = Apr' * rIu(rknot > 0);
                % alpha = (Ap*r)^T/(Ap'*Ap)
                steplength = Ar / ( Apr' * Apr );
                % the right point
                if steplength >= beta
                    steplength = beta;
                    leftb = knotri;
                    left = knoti;
                elseif steplength <= alpha
                    steplength = alpha;
                    rightb = knotri;
                    right = knoti;
                else
                    retcode = [1,0];
                    break;
                end
            end
            if retcode(1) == 1 && retcode(2) == 0
                break;
            end
        end    
        
        if display
            [alpha, minf, knot,retcode] = arraySpiece(A,b,x0,p,tol,maxits);
            [~, normr, ~, g, normKKT, face1, faceN] = kktResidual(A, b, x0, [], 1);
            xt = x0 + steplength * p; xt(xt<0) =0;
            [~, normrN, ~, g0, normKKTN, faceN1, faceN2] = kktResidual(A, b, xt, [], 1);
            pg = g0' * -g; 
            fprintf('simple(steepdown,%d): normr(%g,%g,%g,%g),pg(%g),activeX(%d,%d),activeB(%d,%d) alpha(%g,%g)\n'...
                ,i,normr,normrN,normKKT,normKKTN,pg,face1,faceN1,face2,faceN2,steplength,alpha);
            knoty = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), [testalphax,testalphab]);
            
%             xa = [0 : steplength / 100 : steplength * 1.2];
%             ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
%             pxy={};
%             pxy(1).X = knot;
%             pxy(1).Y = knoty;
%             pxy(1).X = xa;
%             pxy(1).Y = ya;
%             figure
%             hold on
%             p1 = arrayfun(@(a) plot(a.X,a.Y),pxy,'UniformOutput',false);
%             p1(1).Marker = 'o';
%             p1(2).Marker = '+';
%             hold off
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
        switch(newton)
            case 'Lsqr'
                [u,flag,relres,fs,resvec,lsvec,out] = lsqrmx( AI ,bI,lsqrTol,maxIter,[],[],zeros(size(AI,2),1),A(:,RR),b,x0(RR),I,3);
               if norm(u) > 1e-10
%                     x0(:) = 0;
                    x0(RR) = x0(RR) + u;
               end
            case 'ml'
                xz(:) = 0;
                xz(RR) = (AT(RR, I) * A(I,RR))\(AT(RR, I) * b(I));
                r = b - A(:,RR) * xz(RR);
                rb = r( ~I );
                r(r < 0) = 0;
                Ar = AT * r;
                fs = 0;
                if all(xz(RR) > -tol) &  all( rb < tol  ) & all(Ar < tol)
                    fs = 1;
                    x0 = xz;
                end
        end

        if display
            u = xz - x0;
            [rpk, normrN, minx, g, normKKTN, face1N, face2N] = kktResidual(A, b, x0, [], 1);
            % record the value of objection function
            resvec(iter +1) = normrN;
            % record the value of the gradient function
            arvec(iter +1) = normKKTN;
            face1vec(iter +1) = face1N;
            face2vec(iter +1) = face2N;
            pg = g' * (xz-x0);
            [minx,loc]=min(x0);
            fprintf('newton(%d end %d): norm(%g,%g),normKKT(%g,%g),minx(%g,%d),gp(%g,%g),xa(%d,%d),ba(%d,%d)\n',...
                (iter-1)/2,fs,normr,normrN,normKKT,normKKTN,minx,loc,pg,norm(u),face1,face1N,face2,face2N);
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

