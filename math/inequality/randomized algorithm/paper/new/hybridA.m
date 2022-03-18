function [xk,flag,iter,error_k,indexsm] = hybridA(A,b,x0,maxit,nf,type,tol,exactx,debug)

[m,n] = size(A);

%[rpk, r0, normr, normAr] = residual(A,b,x0);
rpk=b-A*x0;
r=rpk;
r(r<0)=0;
iter = 0;

norm_r = norm(r);
if isempty(exactx)
    error_k = [norm(A'*r);];
    iter_k =[0];
else
    error_k = [norm(x-exactx)];
    iter_k =[x];
end

indexsm = 0;
% flag 0-4 return lsqr flag
flag = 5;
[Q,R]=qr(A);
Qn=Q(:,1:n);
% [xfA,rpk] = fmnf(A,b,x0,n,Q,R,rpk,nf);
%while normAr > tol * normA * normr && normr > tol
%while normAr > tol  && normr > tol
for i=1:maxit
    iter = iter + 1;
    [xfA,rpk] = simple(A,b,x0,n,Q,R,rpk,nf);
    isSub = strategies(A,b,Qn,type,iter,nf,rpk,xfA);
    if isSub
        xf = xfA(:, end);
        [xs,rpk,len,flag] = sm(A, b, n, rpk, xf);
        indexsm = indexsm + 1;
        x0 = xs;
    else
        x0 = xfA(:, end);
    end
    if  flag == 0
       fprintf("error flag");
        break;
    end
    if flag < 5 || flag > 6 
%          if flag < 12
%              disp(['flag:',num2str(flag)]);
%          end
         flag = 5;
    end
    r=rpk;
    r(r<0)=0;
        norm_rn = norm(r);
        if norm(A'*r)<tol || norm_rn < 1e-6
            break;
        end
        norm_r = norm_rn;
    end
    % 主要记录迭代过程中的值 用来调试
    if debug
            if ~isempty(exactx)
                e = norm(x-exactx);
                iter_k =[iter_k x];
            else
                Ar = A'*r;
                e = norm(Ar);
                iter_k =[iter_k i];
            end
            error_k = [error_k,e];
            index_k = [index_k,pickedj];
    end
    xk = x0;
end
