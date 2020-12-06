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
[normr,g,normKKT,Ax] = dffunc(A,b,x0);
resvec = zeros(1,maxit);
% the normal gradient
arvec = zeros(1,maxit);
% face b-Ax
faceXvec = zeros(1,maxit);
resvec(1) = normr;
arvec(1) = normKKT;

while norm( x0 .* g, inf) > tol || min( g )< -tol
    ADA = AT(:,I) * Ax(I);
    Ix = x > 1e-10;
    ADA = ;
    %             ADA = A(I,:)' * Ax(I);
    d = x0./(ADA + det);
    fk=fq(A,b,x0);
    AA=det2F(x0,A,b,M)+diag(ones(size(A',1),1)).*delt;
    % 计算方向
    p=-1*AA\d1;
    % 计算步长
    k=k+1;
    
end
toc(t);
x0(x0<0)=0;
x=x0;
fk=fq(A,b,x);
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
