% 网上下载的一个wolf条件实现算法
function [alpha, retcode] = wolfe(A,b,xk, dk, range, maxit)
% alphaMin=t*min(-1.*xk./dk);
rho = 0.25; sigma = 0.75;
[~, normr, ~, Ar, ~ , ~, ~] = kktResidual(A, b, xk , [], []);
qp = rho * Ar' * dk;
alpha = range;
loopcount = 0;
while func(A , b, xk, dk ,alpha) > normr +  alpha * qp   
    loopcount = loopcount + 1;    
    alpha=alpha/2;
    if loopcount > maxit
        alpha = range;
        retcode = [2,loopcount];
        return;
    end
end
loopcount = 0;
while func(A , b, xk, dk ,alpha) > normr +  alpha * qp
        alpha=alpha/2;
end
beta = 2 * alpha;
if ~(fdetq(A,b,xk+alpha*dk)'*dk >= sigma*fdetq(A,b,xk)'*dk)
%         aa = alpha;
%         alpha = min([2*alpha, (bb+alpha)/2]);
        alpha=alpha/2;
        continue;
    end
%     if alpha<0 || alpha>alphaMin
%         index=index+1;
%         alpha=alphaMin-index*detAlpha;
%         continue;
%     end
    break;


newxk = xk+alpha*dk;
newfdetqk = fdetq(A,b,newxk);
newfk = fq(A,b,newxk);
end 
function fvalue = func(A,b,x0,p,alpha)
r = b - A * (x0 + alpha * p);
r( r < 0 ) = 0;
fvalue = 0.5 * (r' * r);
end
