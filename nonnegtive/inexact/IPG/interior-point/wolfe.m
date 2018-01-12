function [alpha, newxk, newfk, newfdetqk] = wolfe(A,b,xk, dk, t)
alphaMin=t*min(-1.*xk./dk);
rho = 0.25; sigma = 0.75;detAlpha=0.001;
alpha = alphaMin; a = 0; b = Inf; 
index=1;
while (1)
    if ~(fq(A,b,xk+alpha*dk)<=fq(A,b,xk)+rho*alpha*fdetq(A,b,xk)'*dk)
        b = alpha;
        alpha = (alpha+a)/2;
        continue;
    end
    if ~(fdetq(A,b,xk+alpha*dk)'*dk >= sigma*fdetq(A,b,xk)'*dk)
        a = alpha;
        alpha = min([2*alpha, (b+alpha)/2]);
        continue;
    end
    if alpha<0 || alpha>alphaMin
        index=index+1;
        alpha=alphaMin-index*detAlpha;
        continue;
    end
    break;
end
newxk = xk+alpha*dk;
newfdetqk = fdetq(A,b,newxk);
newfk = fq(A,b,newxk);