% 罚函数方法
% （海森矩阵+防止错误项）* 方向 = 梯度
%  搜索步长
function [xk, resvec, arvec, faceXvec, tf] = GNP(A,b,x0,M,delt,tol,maxit)
t=clock;
loopcount = 1;
%[rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
% the residual vector
AT = A';
[m,n] =size(A);
MD = speye(n)*M;
Ddelt = speye(n)*delt;
[normr,g,normKKT,Ax] = dffunc(A,b,x0);
resvec = zeros(1,maxit);
% the normal gradient
arvec = zeros(1,maxit);
% face b-Ax
faceXvec = zeros(1,maxit);
resvec(1) = normr;
arvec(1) = normKKT;

while norm( x0 .* g, inf) > tol || min( g )< -tol
    I = b > Ax;
    ADA = AT(:,I) * A(I,:);
    Ix = x0 > 1e-10;
    MD0 = MD;
    MD0(Ix,Ix) = 0;
    ADA = ADA + 2 * MD0 + Ddelt;
    % ADA = A(I,:)' * Ax(I);
    % 计算方向
 %   p = lsqminnorm(ADA, -g);
    p = -ADA \ g;
   [alpha, knot, retcode] = arraySpiece(A,b,x0,p,1e-5,30);
    % 计算步长
    loopcount=loopcount+1;
    x0 = x0 + alpha * p;
    x0(x0<0) = 0;
    [normr,g,normKKT,Ax] = dtffunc(A,AT,b,x0);
    %    [normr,g,normKKT,Ax] = dffunc(A,b,x0);
%     [rpk, normr, ~, g, normKKT, faceX, ~] = kktResidual(A, b, x0 , [], 1);
    % record the value of objection function
    resvec(loopcount) = normr;
    % record the value of the gradient function
    arvec(loopcount) = normKKT;
end
tf = etime(clock,t);
x0(x0<0)=0;
xk = x0;
end

function [fvalue,gvalue,normKKT,Ax] = dtffunc(A,AT,b,x0)
Ax = A * x0;
r = b - Ax;
r( r < 0 ) = 0;
gvalue =  -AT * r;
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
