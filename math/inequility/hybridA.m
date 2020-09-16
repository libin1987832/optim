function [xk,flag,relres,iter,resvec,arvec,itersm,tf] = hybridA(A,b,x0,maxit,nf,type)
t=clock;

% stop criterion
tol = 1e-15;
[m,n] = size(A);

% if n/m < 0.81
%     tol = 1e-12;
% else
%     tol = 1e-10;
% end

normA = norm(A,2);
rpk = b-A * x0;
r0 = rpk;
r0(r0<0) = 0;
normAr = norm(A' * r0);
normr = norm(r0);
iter = 0;
% residual vector
resvec = zeros(1,maxit + 1);
% the normal gradient 
arvec = zeros(1,maxit + 1);
% subspace minization
itersm = zeros(1,maxit + 1);
resvec(1) = normr;
indexsm = 0;
% flag 0-4 return lsqr flag
flag = 5;
[Q,R]=qr(A);
Qn=Q(:,1:n);
% [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
%while normAr > tol * normA * normr && normr > tol
while normAr > tol  && normr > tol
    iter = iter + 1;
    [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
    isSub = strategies(A,b,Qn,iter*nf,type,rpk,xfA);
    if isSub
        xf = xfA(:, end);
%        for j=1:5
        [xs,rpk,len,flag] = sm(A, b, n, rpk, xf);
 %        xf=xs;
 %       end
        indexsm = indexsm + 1;
        itersm(iter + 1) = len;
        x0 = xs;
    else
 %       [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
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
        if flag < 12
            disp(['flag:',num2str(flag)]);
        end
        flag = 5;
    end
end
xk = x0;
relres = normr;
resvec = resvec(1:iter + 1);
itersm = itersm(1:iter + 1);
arvec = arvec(1:iter + 1);
tf = etime(clock,t);