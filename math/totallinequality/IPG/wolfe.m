% 网上下载的一个wolf条件实现算法
function [alpha, newxk, newfk, newfdetqk] = wolfe(A,b,xk, dk, t)
% alphaMin=t*min(-1.*xk./dk);
alphaMin=1;
n=size(xk,1);
for i=1:n
    if dk(i)<0
        ta=-1*xk(i)/dk(i);
        if ta<alphaMin
            alphaMin=ta;
        end
    end
end
alphaMin=t*alphaMin;
rho = 0.25; sigma = 0.75;detAlpha=0.001;
alpha = alphaMin; aa = 0; bb = Inf; 
index=1;
while (1)
    if ~(fq(A,b,xk+alpha*dk)<=fq(A,b,xk)+rho*alpha*fdetq(A,b,xk)'*dk)
        
        alpha=alpha/2;
        continue;
    end
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
end
newxk = xk+alpha*dk;
newfdetqk = fdetq(A,b,newxk);
newfk = fq(A,b,newxk);
