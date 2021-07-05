function [xs,rpk,len,flag]=sm(A,b,n,rpk,x0)
% lsqr tollera
tol = 1e-15;
lsqrTol = 1e-15;
maxIter = 5;

AA = (rpk>tol);
AI = A(AA,:);
bI = rpk(AA);

[u,flag,relres,iter,resvec,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA);
if 0
    xs = x0 + u;
    len = iter;
    rpk = b - A * xs;
else
%   u = AI\bI;
     u=zeros(n,1);
 [U,S,V]=svds(AI);
  [m,n] =size(AI);
    minmn = min(size(S));
 svdj=minmn;
 for j=1:minmn
     if S(j,j)<1e-20
         svdj=j;
         break;
     end
 end
 B=U(:,1:svdj)'*rpk(AA);
 for j=1:svdj
     B(j)=B(j)/S(j,j);
 end
u=V(:,1:svdj)*B;
         aa = spiecewise(A,b,u,x0);
    %aa = spiecewise(A,b,u,x0);
    xs = x0 + aa * u;
    rpk = b - A * xs;
    x0=xs;
    for i=1:6
        I = find(rpk>=tol);
        AI = A(I,:);
        hk = AI \ rpk(I);
        aa = spiecewise(A,b,hk,x0);
        xs = x0 + aa * hk;
        x0 = xs;
        rpk = b - A * x0;
    end
    len = aa;
    flag = flag + 6;
end
%rpk = b - A * xs;
%semilogy(1:iter,resvec(1:iter),'b.');
% if flag == 0 && ~out
%     flag = 0;
% else
%     flag = flag+10;
% end


% normr0 = norm(bI);
% rk=rpk;
% rk(rk<0)=0;
% normrk = norm(rk);
%
% fprintf('sm:subspace:%d iter:%d len:%g r0:%g rk:%g u:%g \n',sum(AA),iter,len,normr0,normrk,norm(u));