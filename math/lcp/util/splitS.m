function [x,res] = splitS(Q,d,s,x0,iter)
addpath('./test')
[m,n]=size(Q);
Di=diag(Q);
Di(Di<0)=0;
% Di=1./Di;
tol=1e-12;
% x0=zeros(n,1);

for i=1:iter
%     if issparse(Q)
%         x=sparse(n,1);
%     else
%         x=zeros(n,1);
%     end
        x=x0;
        for j=1:n
            old_xi = x(j);
            ri     = d(j) + Q(j,:)*x;
            Aii    = Q(j,j);           
            x(j) = max( 0, old_xi - (ri / Aii) );
        end
%     for j=1:n
% %         sd=s*Di(j)*(d(j)+Q(j,:)*[xk(1:(j-1));x0(j:n)]);
%         sd=s*(d(j)+Q(j,:)*[xk(1:(j-1));x0(j:n)])/Di(j);
%         if x0(j)>sd
%             xk(j)=x0(j)-sd;
%         end
%     end
%     xk(xk<0)=0;
    x0=x;
    res=test_valid(Q,d,x);
    if res<tol
        break;
    end
end