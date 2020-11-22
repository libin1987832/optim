addpath('../algorithmInequality')
addpath('../dataInequality')
addpath('../resultInequality')
clear
clc
m = 1000;
n = 500;
[A,b,x0]=randInequality(m,n,2,-1);

iter = 20;
res = zeros(6,iter);
rk0 = b - A * x0;
tic
%D = zeros(n,1);
D = A.*A;
D = sum(D,1);
rk = rk0;
xk = x0;
for i = 1:iter
    [rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
    res(1,i) = nRkN;
    [xk,rk] = FixedGS(xk,A,b,D,rk,3);
end
toc

tic
rk = rk0;
xk = x0;
for i = 1:iter
    [rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
    res(2,i) = nRkN;
    [xk,rk]=FMD(xk,A,b,rk,0);
end
toc

tic
rk = rk0;
xk = x0;
for i = 1:iter
    [rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
    res(3,i) = nRkN;
    [xk,rk]=FMD(xk,A,b,rk,1);
end
toc

tic
rk = rk0;
xk = x0;
[Q,R] = qr(A);
for i = 1:iter
    [rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
    res(4,i) = nRkN;
    [xk,rk]=FMQR(xk,Q,R,A,b,rk);
end
toc

tic
xk = x0;
rk = rk0;
[Q,R] = qr(A);
for i = 1:iter
    [rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
    res(5,i) = nRkN;
    [xk,rk]=Lei(xk,A,b,3,rk);
end
toc

tic
xk = x0;
rk = rk0;
ATA = A' * A;
D = diag(diag(ATA));
L = tril(ATA) - D;
U = triu(ATA) - D;
iDL = inv(D+L);
for i = 1:iter
    [rk, rkN, nRkN, nGARk] = residual(A,b,xk,rk);
    res(6,i) = nRkN;
    zk = -rk;
    zk(zk<0) = 0;
    ATb = A' * ( b + zk );
    xk = -iDL * ( U * xk -  ATb);
    rk = b - A * xk;
end
toc


type=['r','c','k','g','w','b'];
typet=['+','o','v','s','.','+'];
plotSemilogy(res,type,typet);