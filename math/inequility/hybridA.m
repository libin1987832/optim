function [xk,flag,relres,iter,resvec,itersm,tf] = hybridA(A,b,x0,maxit,nf,type)
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
resvec = zeros(1,maxit+1);
itersm = false(1,maxit+1);
resvec(1) = normr;
flag = 0;
[Q,R]=qr(A);
Qn=Q(:,1:n);
while normAr > tol * normA * normr;
    iter = iter + 1;
    [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
    isSub = strategies(A,b,Qn,iter*nf,type,rpk,xfA);
    itersm(iter + 1) = isSub;
    if isSub
        xf = xfA(:, end);
        xs = sm(A, b, n, rpk, xf);
        x0 = xs;
    else
        x0 = xfA(:, end);
    end
    rpk = b - A * x0;
    r = rpk;
    r(r<0) = 0;
    normAr = norm(A' * r);
    resvec(iter + 1)=normAr;
    normr = norm( r );
    if iter > maxit
        flag = 1;
        break;
    end
end
xk = x0;
relres = normr;
resvec = resvec(1:iter + 1);
itersm = itersm(1:iter + 1);
tf = etime(clock,t);