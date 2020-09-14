function [x,flag,relres,iter,resvec] = hybridA(A,b,x0,maxit,type)
normA = norm(A,2);
rpk = b-A * x0;
r0 = rpk;
r0(r0<0) = 0;
normAr = norm(A' * r0);
normr = norm(r0);
iter = 0;
resvec = zeros(maxit,1);
resvec(1) = normr;
flag = 0;
[Q,~]=qr(A);
Qn=Q(:,1:n);
while normAr > tol * normA * normr;
    iter = iter + 1;
    xf = fmnf(A,b,Qn,x0,rpk,nf);
    isSub = strategies(A,b,Qn,type);
    if isSub
        xs = sm(A,b,xf);
        x0 = xs;
    else
        x0 = xf;
    end
    rpk = b - A * x0;
    r = rpk;
    r(r<0) = 0;
    normAr = norm(A' * r);
    resvec(iter+1)=normAr;
    normr = norm( r );
    if iter > maxit
        flag = 1;
        break;
    end
end
relres = normr;
resvec = resvec(1:iter+1);