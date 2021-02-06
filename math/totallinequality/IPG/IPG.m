%[xk,resvec,arvec,face1vec,face2vec,tf]
function [xk, resvec, arvec, faceXvec, tf] = IPG(A,b,x0,tol,det,tou,maxit,dtype)
display = 1;
%x
t=clock;
loopcount = 1;
%[rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
% the residual vector
AT = A';
[m,n] = size(A);
[normr,g,normKKT,Ax] = dffunc(A,b,x0);
resvec = zeros(1,maxit);
% the normal gradient
arvec = zeros(1,maxit);
% face b-Ax
faceXvec = zeros(1,maxit);
resvec(1) = normr;
arvec(1) = normKKT;
faceXvec(1) = 0;
IF = false(m,1);
IP = b - Ax > 1e-15;
ADA = AT(:,IP) * Ax(IP);
% dtype = 'NT';
while norm( x0 .* g, inf) > tol || min( g )< -tol
    if loopcount > maxit
        break;
    end
    loopcount = loopcount +1;
    
    switch dtype
        case 'NT'
                         p = -(A(IP,:)'*A(IP,:))\g;
            %p = lsqminnorm(A(IP,:)'*A(IP,:), -g);
            %  cos2 =  np'*p / (norm(np)*norm(p));
        case 'ST'
            p = -g;
        otherwise
            %           ADA = A(I,:)' * Ax(I);
            d = x0./(ADA + det);
            p = - d .* g;
    end
    Ap = A * p;
    alphaAll = - x0./p;
    alphak = min(alphaAll(alphaAll>0));
    if isempty(alphak)
        alphak = 1;
    end
    %tou = toustrategy(tou,loopcount,x0);
    % tou = tou - 1/100;
    touk = tou + rand()*(1-tou);
    talphak = touk * alphak;
    % assert descent direction
    
     [alpha, knot, retcode] = arraySpiece(A,b,x0,p,1e-10,30);
      x00 =x0;    g00=g; Ax00=Ax;
%     [alpha, ~, retcode] = wolfe(A, b, x0, p, talphak,Ax,Ap,normr,g,100);
    
    x0 = x0 + alpha * p;
    Ax = Ax + alpha * Ap;
    r = b - Ax;
    IP =  r > 1e-15;
    ADA = AT(:,IP) * Ax(IP);
    if any(xor(IF,IP))
        ADb = AT(:,IP) * b(IP);
    end

    g = ADA - ADb;
%     r=b-A*x0;
%     r(r<0)=0;
%     gr = -A'*r;
    IF = IP;
    normr = 0.5*r(IP)'*r(IP);
    %    [normr,g,normKKT,Ax] = dtffunc(A,AT,b,x0);
    
    if display
        [normrd,g1,normKKT11,Ax1] = dtffunc(A,AT,b,x00);
        [normr0,g10,normKKT10,Ax10] = dtffunc(A,AT,b,x0);

        pg = g1'*p;% assert descent direction
        [minx,loc]=min(x0);
        fprintf('IPG(%d end): normr(%g),normrk(%g),kkt(%g),gp(%g),alpha(%g),talpha(%g),minx(%g,%d)\n'...
            , loopcount, normrd, normr0,normKKT11,pg,alpha,alphak,minx,loc);
       if normr0 > normrd
        dd
        end
            resvec(loopcount) = normr;
    % record the value of the gradient function
    arvec(loopcount) = normKKT11;
    faceXvec(loopcount) = alphak;
    end

    
end
xk=x0;
tf = etime(clock,t);
end

function [fvalue,gvalue,normKKT,Ax] = dtffunc(A,AT,b,x0)
Ax = A * x0;
r = b - Ax;
r(r<0) = 0;
gvalue =  - AT * r;
normKKT = norm( x0 .* gvalue, inf);
fvalue = 0.5 * (r' * r);
end
function [fvalue,gvalue,normKKT,Ax] = dffunc(A,b,x0)
Ax = A * x0;
r = b - Ax;
r( r < 0 ) = 0;
gvalue =  -A' * r;
normKKT = norm( x0 .* gvalue, inf);
fvalue = 0.5 * (r' * r);
end

function fvalue = func(A,b,x0)
Ap = A * x0;
r = b - Ap;
r( r < 0 ) = 0;
fvalue = 0.5 * (r' * r);
end

function fvalue = dfunc(A,b,x0)
Ap = A * x0;
r = b - Ap ;
r( r < 0 ) = 0;
fvalue =  -A' * r;
end

%% test example
% A = [1 3;2 4;-5 -6]; b = [5;6;-3]; x0 = [1;1];
% -g = [39;45]; ADA = [26 33;33 45]; ADAx = [59;78];
% p = [39/59;45/78] pnt = [30/9; -13/9]
% alpha = +inf alpha = 9/13  alpha =1 f = 8+39*5/59+180/78=13.6128^2

