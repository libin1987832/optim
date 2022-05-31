function [x,k,R]=GSO_improve(A,b,n,RRE,x,b_A)
%b_A=b-(A'\(A'*b));
RRE=RRE*norm(b,2)^(2);
N=zeros(n,1);%用于存储矩阵列范数
G=zeros(n,1);
R=[];
for i=1:1:n-1
    N(i)=norm(A(:,i),2)^(2);
    G(i)=A(:,i)'*A(:,i+1);
end
N(n)=norm(A(:,n),2)^(2);
G(n)=A(:,n)'*A(:,1);
GN=-G./N;
g=GN.*G+circshift(N,-1);
%第一次迭代
 r=b-A*x;
 R=r;
 alpha=A(:,1)'*r/N(1);
 x(1)=x(1)+alpha;
 r=r-alpha*A(:,1);
 i_1=1;
for k=1:1:1000
     i_2=rem(k,n)+1;
     if g(i_1)>0.0001
         alpha=A(:,i_2)'*r/g(i_1);
         beta=GN(i_1)*alpha;
         x(i_2)=x(i_2)+alpha;
         x(i_1)=x(i_1)+beta;
         r=r-alpha*A(:,i_2)-beta*A(:,i_1);
         R=[R,r];
     end
     if  norm(r-b_A)^2<RRE
         break;
     else
         i_1=i_2;
    end
end
k=k+1;
end
