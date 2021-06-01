function [xk,rk] = FMGS2(x0,A,b,D,rk,niter)
[m,n]=size(A);
xk=x0;
zk=-rk;
zk(zk<0)=0;
bzk=b+zk;
% compute min increase
for j=1:niter
    for i=1:n
      r=bzk-A*xk;
       xk(i) = xk(i) + (sum(A(:,i).*r))/D(i);
%xk(i) = xk(i) + (sum(A(:,i).*(bzk-A*xk)))/D(i);
    end
end
rk=b-A*xk;
% function x = FMGS2(A, b, x, iters)
%     for i = 1:iters
%         for j = 1:size(A,1)
%             x(j) = (b(j) - sum(A(j,:).*x) + A(j,j)*x(j))/A(j,j);
%         end
%     end
% end