
% compuation lsqrm: A*x(n)   lsqminnorm: A*x(n) A*xs spiecewise A*x(n)
function [xs,rpk,len,flag]=sm(A,b,n,rpk,x0)
% lsqr tollera
tol = -1e-15;
lsqrTol = 1e-15;
maxIter = 5;
% AA indice of active set
I = (rpk>tol);
% subspace
AI = A(I,:);
bI = rpk(I);

% [u,flag,relres,iter,resvec,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA);
% if ~out
%     xs = x0 + u;
%     len = 3;
%     rpk = b - A * xs;
% else
%    u = AI\bI;

%      u = lsqminnorm(AI,bI);

%     aa = spiecewise(A,b,u,x0);
%     xs = x0 + aa * u;
%     rpk = b - A * xs;
%     x0=xs;

   for i=1:1


 %       hk = lsqminnorm(AI,rpk(I));
 %      hk = AI \ rpk(I);
 hk=zeros(n,1);
 [U,S,V]=svd(AI);
 minmn=min(size(S));
 svdj=minmn;
 for j=1:minmn
     if S(j,j)<1e-10
         svdj=j;
         break;
     end
 end
 B=U(:,1:svdj)'*rpk(I);
 for j=1:svdj
     hk(j)=B(j)/S(j,j);
 end
 hk=V(:,1:svdj)*hk;
aa = spiecewise(A,b,hk,x0);
        xs = x0 + aa * hk;
        x0 = xs;
        rpk = b - A * x0;
        r=rpk;
        r(r<0)=0;
        nr=A'*r;
        
  % end
    len = aa;
    flag = 1+ 6;
   end
% end
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