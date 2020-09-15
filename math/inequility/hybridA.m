function [xk,flag,relres,iter,resvec,arvec,itersm,tf] = hybridA(A,b,x0,maxit,nf,type)
t=clock;
tol = eps;
[m,n] = size(A);
normA = norm(A,2);
rpk = b-A * x0;
r0 = rpk;
r0(r0<0) = 0;
normAr = norm(A' * r0);
normr = norm(r0);
iter = 0;
resvec = zeros(1,maxit + 1);
arvec = zeros(1,maxit + 1);
itersm = zeros(1,maxit + 1);
resvec(1) = normr;
indexsm = 0;
flag = 5;
[Q,R]=qr(A);
Qn=Q(:,1:n);
while normAr > tol * normA * normr;
    iter = iter + 1;
    [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
    isSub = strategies(A,b,Qn,iter*nf,type,rpk,xfA);
    if isSub
        xf = xfA(:, end);
        [xs,rpk,len,flag] = sm(A, b, n, rpk, xf);
        indexsm = indexsm + 1;
        itersm(iter + 1) = len;
        x0 = xs;
    else
        itersm(iter + 1) = -1;
        x0 = xfA(:, end);
    end
    r = rpk;
    r(r<0) = 0;
    normAr = norm(A' * r);
    normr = norm( r );
    resvec(iter + 1) = normr;
    arvec(iter + 1) = normAr;
    if iter > maxit || flag == 0
        break;
    end
    if flag < 5 || flag > 6
        disp(['flag:' flag]);
    end
end
xk = x0;
relres = normr;
resvec = resvec(1:iter + 1);
itersm = itersm(1:iter + 1);
arvec = arvec(1:iter + 1);
tf = etime(clock,t);