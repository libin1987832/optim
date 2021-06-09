% flag
% relres
% iter
% resvec
% arvec
% itersm
% tf
function [xk,flag,relres,iter,resvec,arvec,itersm,tf] = hybridA(A,b,x0,steplength,maxit,nf,type)
t=clock;

% stop criterion
tol = 1e-13;
[m,n] = size(A);


%normA = norm(A,2);
% rpk = b-A * x0;
% r0 = rpk;
% r0(r0<0) = 0;
% normAr = norm(A' * r0);
% normr = norm(r0);

[rpk, r0, normr, normAr] = residual(A,b,x0);

iter = 0;
% residual vector
resvec = zeros(1,maxit + 1);
% the normal gradient 
arvec = zeros(1,maxit + 1);
% subspace minization
itersm = zeros(1,maxit + 1);
resvec(1) = normr;
arvec(1) = normAr;
indexsm = 0;
% flag 0-4 return lsqr flag
flag = 5;
[Q,R]=qr(A);
Qn=Q(:,1:n);
% [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
%while normAr > tol * normA * normr && normr > tol
while normAr > tol  && normr > tol
    iter = iter + 1;
%    [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf,type);
    [xfA,rpk] = simple(A,b,x0,n,Q,R,rpk,nf,type);
%     debug = 0;
%     if debug
%        fprintf("active:%d\n",sum(sign(rpk)>0));
%     end
    isSub = strategies(A,b,Qn,type,iter,nf,rpk,xfA);
    if isSub
        xf = xfA(:, end);
        iter
        [xs,rpk,len,flag] = sm(A, b, n, rpk, xf);
        indexsm = indexsm + 1;
        % record step length and statistic to number of sm
        itersm(iter + 1) = len;
        x0 = xs;
    else
 %       [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
        itersm(iter + 1) = -1;
        x0 = xfA(:, end);
    end
    [rpk, r, normr, normAr] = residual(A,b,x0,rpk);
    % record the value of objection function
    resvec(iter + 1) = normr;
    % record the value of the gradient function
    arvec(iter + 1) = normAr;
    if iter > maxit || flag == 0
        break;
    end
    if flag < 5 || flag > 6 
%         if flag < 12
%             disp(['flag:',num2str(flag)]);
%         end
         flag = 5;
    end
end
xk = x0;
relres = normr;
resvec = resvec(1:iter + 1);
itersm = itersm(1:iter + 1);
arvec = arvec(1:iter + 1);
tf = etime(clock,t);