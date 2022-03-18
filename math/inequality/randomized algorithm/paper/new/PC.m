function [xk,zk,iter,error_k,index_k]=PC(x1,z1,A,b,maxit,tol,exactx,debug)
[m, n] = size(A);
% x1 = x0;
% z1=z0
Ax1 = A*x1;
r=b-Ax1;
rp=r;
rp(rp<0)=0;
norm_r=norm(rp);
iter = 0;
if isempty(exactx)
    error_k = [norm(A'*rp)];
    iter_k =[0];
else
    error_k = [norm(x-exactx)];
    iter_k =[x];
end
index_k=[0];
for i = 1:maxit
    z=-r;
    e1=A'*(z-z1);
  %  e1=A'*(A*x1-z1-b);
  %  z=A*x1-b;
    z(z<0)=0;
    e2=z1-z;
    sum=e1'*e1+e2'*e2;
    ae=A*e1-e2;
    p=(sum)/(sum+ae'*ae);
    xk=x1-p*e1;
    zk=z1-p*e2;
    x1=xk;
    r=b-A*x1;
    z1=zk;
    rp=r;
    rp(rp<0)=0;
    norm_rn=norm(rp);
    iter=iter+1;
    
    if abs(norm_rn-norm_r)<tol || norm_rn < 1e-6
        break;
    end
    norm_r=norm_rn;
    % 主要记录迭代过程中的值 用来调试
    if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            iter_k =[iter_k x];
        else
            Ar = A'*rp;
            e = norm(Ar);
            iter_k =[iter_k i];
        end
        error_k = [error_k,e];
        index_k = [index_k,i];
    end
    
end
