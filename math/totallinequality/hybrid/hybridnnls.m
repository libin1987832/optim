function [xk,resvec,arvec,tf] = hybridnnls(A,b,x0,alpha,nf,maxit)
t=clock;

% stop criterion
tol = 1e-13;
[m,n] = size(A);
lsqrTol = 1e-13;
maxIter = 20;
[rpk, normr, xmin, Ar, normKKT] = kktResidual(A, b, x0 , 1);
iter = 0;
% the residual vector
resvec = zeros(1,maxit + 1);
% the normal gradient 
arvec = zeros(1,maxit + 1);
% subspace minization
itersm = zeros(1,maxit + 1);
resvec(1) = xmin;
arvec(1) = normKKT;
indexsm = 0;
% flag 0-4 return lsqr flag
flag = 5;
xfA = zeros(n,nf);
xfA(:,1) = x0;
% [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
%while normAr > tol * normA * normr && normr > tol
while normKKT > tol 
    iter = iter + 1;
    for i = 2:nf+1
    % AA indice of active set
    AA = (rpk>tol);
    % subspace
    AI = A(AA,:);
    bI = rpk(AA);
    
    [u,flag,relres,iter,resvec,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA,2);
    xls = x0 + u;
    [rpk, normr, xmin, Ar, normKKT] = kktResidual(A, b, xls , type)
    rpk = b - A * xls;
    [x,rpk] = projectfixed(A,b,xls,rpk,alpha);
    [rpk, normr, xmin, Ar, normKKT] = kktResidual(A, b, x , type)
    xfA(:,i) = x;
    end
    isSub = strategy(A,b,steplengthOrk,2,iter,nf,rpk,xfA);
    if isSub
        xf = xfA(:, end);
        AA = (rpk>tol);
        % subspace
        AI = A(AA,:);
        bI = rpk(AA);
        [u,flag,relres,iter,resvec,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA,3);
        indexsm = indexsm + 1;
        % record step length and statistic to number of sm
        itersm(iter + 1) = len;
        x0 = xs;
    else
 %       [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
        itersm(iter + 1) = -1;
        x0 = xfA(:, end);
    end
%     r = rpk;
%     r(r<0) = 0;
%     normAr = norm(A' * r);
%     normr = norm( r );
    [rpk, r, normr, normAr] = residual(A,b,x0,rpk);
    % record the value of objection function
    resvec(iter + 1) = normr;
    % record the value of the gradient function
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